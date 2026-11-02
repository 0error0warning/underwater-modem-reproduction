function report = validate_packet_link_v4()
%VALIDATE_PACKET_LINK_V4 Framing, CRC, blind waveform and stop-wait boundaries.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'matlab'));
out=fullfile(root,'results','v4','packet_link'); if ~exist(out,'dir'), mkdir(out); end
cfg=uaclink.default_config(); rows={};
rows{end+1}=test('crc_standard_123456789',double(uaclink.crc16(uint8('123456789')))==hex2dec('29B1'),'CCITT-FALSE check=0x29B1');
rows{end+1}=test('crc_empty',double(uaclink.crc16(zeros(1,0,'uint8')))==65535,'Empty input keeps initFFFF');
message='你好，水声 Modem!'; payload=unicode2native(message,'UTF-8');
[wire,bits]=uaclink.encode_frame(payload,257,1,cfg); parsed=uaclink.decode_frame(wire,cfg);
rows{end+1}=test('utf8_wire_roundtrip',parsed.valid && parsed.sequence==257 && ...
    strcmp(native2unicode(parsed.payload,'UTF-8'),message) && isequal(wire,uaclink.bits_to_bytes(bits)),'Header plus UTF8 payload roundtrip');
acceptedCorrupt=0;
for k=1:numel(bits)
    bad=bits; bad(k)=1-bad(k); p=uaclink.decode_frame(uaclink.bits_to_bytes(bad),cfg);
    acceptedCorrupt=acceptedCorrupt+p.valid;
end
rows{end+1}=test('all_single_bit_corruption',acceptedCorrupt==0,sprintf('%d individual bit flips, header/payload/CRC included',numel(bits)));
for len=[0,64]
    bytes=uint8(0:len-1); [w,~]=uaclink.encode_frame(bytes,65535,1,cfg); p=uaclink.decode_frame(w,cfg);
    rows{end+1}=test(sprintf('payload_length_%d',len),p.valid && isequal(p.payload,bytes),'Boundary wire roundtrip'); %#ok<AGROW>
end
rows{end+1}=test('reject_oversize_encoder',throws(@()uaclink.encode_frame(zeros(1,65),0,1,cfg)),'65 bytes exceeds maximum64');
[ackWire,~]=uaclink.encode_frame(zeros(1,0,'uint8'),65535,2,cfg); p=uaclink.decode_frame(ackWire,cfg);
rows{end+1}=test('empty_ack',p.valid && p.type==2 && isempty(p.payload),'CRC protected zero-length ACK');
for k=1:7
    w=wire;
    reasons={'bad_magic','unsupported_version','unknown_type','length_exceeds_limit','ack_with_payload','frame_length_mismatch','crc_mismatch'};
    switch k
        case 1, w(1)=0;
        case 2, w(3)=2;
        case 3, w(4)=3;
        case 4, w(7:8)=uint8([0,65]);
        case 5, w(4)=2;
        case 6, w=w(1:end-1);
        case 7, w(9)=bitxor(w(9),uint8(1));
    end
    p=uaclink.decode_frame(w,cfg);
    rows{end+1}=test(['reject_' reasons{k}],~p.valid && strcmp(p.reason,reasons{k}) && isempty(p.payload),'No delivery on malformed wire frame'); %#ok<AGROW>
end
% PHY tests exercise the bounded header parser without supplying true length.
definitions={message,0,Inf;'',0.5,20;repmat('X',1,64),-0.5,20;'waveform check',1.5,10;'negative drift',-1.5,10};
for k=1:size(definitions,1)
    txCfg=cfg; txCfg.phase_rad=0.23+0.77*k;
    expected=unicode2native(definitions{k,1},'UTF-8');
    [tx,~]=uaclink.transmit(expected,k-1,1,txCfg);
    rx=uaclink.synthetic_channel(tx,cfg,500+k,definitions{k,2},definitions{k,3});
    p=uaclink.receive(rx,cfg);
    rows{end+1}=test(sprintf('blind_waveform_%d',k),p.valid && isequal(p.payload(:),uint8(expected(:))), ...
        sprintf('bytes=%d v=%+.1fm/s sampleSNR=%g dB reason=%s',numel(expected),definitions{k,2},definitions{k,3},p.reason)); %#ok<AGROW>
end
badBits=bits; badBits(8*8+1)=1-badBits(8*8+1);
corruptTx=uac2.make_packet(badBits,cfg);
p=uaclink.receive(uaclink.synthetic_channel(corruptTx,cfg,700,0,Inf),cfg);
rows{end+1}=test('waveform_crc_corruption',~p.valid && strcmp(p.reason,'crc_mismatch'),'Valid preamble/header but one payload wire bit intentionally changed');
badLength=wire; badLength(7:8)=uint8([0,65]); badTx=uac2.make_packet(uaclink.bytes_to_bits(badLength),cfg);
p=uaclink.receive(uaclink.synthetic_channel(badTx,cfg,701,0,Inf),cfg);
rows{end+1}=test('waveform_length_bound',~p.valid && strcmp(p.reason,'length_exceeds_limit'),'Refuses parsed length before allocating/decode of full frame');
[tx,~]=uaclink.transmit(payload,0,1,cfg); layout=uac2.packet_layout(cfg,numel(bits));
partial=tx(1:layout.payload_start_sample-1+(64+8)*round(cfg.symbol_s*cfg.fs));
p=uaclink.receive([zeros(round(.2*cfg.fs),1);partial],cfg);
rows{end+1}=test('waveform_truncated_body',~p.valid && strcmp(p.reason,'truncated_frame'),'Header intact, declared body unavailable');
p=uaclink.receive(tx(1:round(.1*cfg.fs)),cfg);
rows{end+1}=test('waveform_truncated_preamble',~p.valid,'Only100ms received');
p=uaclink.receive(zeros(round(2*cfg.fs),1),cfg);
rows{end+1}=test('silence_rejected',~p.valid,'Two-second silence');
rng(901,'twister'); falseAccept=0;
for k=1:10
    p=uaclink.receive(randn(round(2*cfg.fs),1),cfg); falseAccept=falseAccept+p.valid;
end
rows{end+1}=test('noise_only_rejected',falseAccept==0,'Zero accepted of10 synthetic noise-only windows; finite-sample check only');
% Node semantics tested through real waveform ACK/DATA decoders.
a=uaclink.new_node('A'); b=uaclink.new_node('B');
[a,tx]=uaclink.send_message(a,'one',cfg,0); deadline=a.deadline_s;
[early,noWave]=uaclink.on_timeout(a,cfg,deadline-1e-6);
rows{end+1}=test('no_early_retry',early.attempts==1 && isempty(noWave),'Timeout does not fire before deadline');
[wrongAck,~]=uaclink.transmit(zeros(1,0,'uint8'),1,2,cfg);
[unchanged,~,event]=uaclink.on_waveform(a,padded(wrongAck,cfg),cfg,.1);
rows{end+1}=test('wrong_sequence_ack',strcmp(unchanged.state,'waiting_ack') && strcmp(event.name,'unexpected_ack'),'ACKseq1 cannot clear pendingseq0');
[~,ackBits]=uaclink.encode_frame(zeros(1,0,'uint8'),0,2,cfg); ackBits(end)=1-ackBits(end);
badAck=uac2.make_packet(ackBits,cfg);
[unchanged,~,event]=uaclink.on_waveform(a,padded(badAck,cfg),cfg,.2);
rows{end+1}=test('corrupted_ack_keeps_pending',strcmp(unchanged.state,'waiting_ack') && ...
    ~event.packet_valid && strcmp(event.reason,'crc_mismatch'),'ACK with one damaged CRC bit cannot clear pendingDATA');
[b,ack,event]=uaclink.on_waveform(b,padded(tx,cfg),cfg,1.6);
rows{end+1}=test('first_data_delivery',event.delivered && numel(b.delivered)==1 && ~isempty(ack),'Data emitted one application payload and modulated ACK');
[b,ackDuplicate,event]=uaclink.on_waveform(b,padded(tx,cfg),cfg,1.7);
rows{end+1}=test('duplicate_delivery_suppressed',~event.delivered && numel(b.delivered)==1 && ~isempty(ackDuplicate),'Lost-ACK retry receives ACK but no second application delivery');
[conflictTx,~]=uaclink.transmit(uint8('different'),0,1,cfg);
[b,noAck,event]=uaclink.on_waveform(b,padded(conflictTx,cfg),cfg,1.8);
rows{end+1}=test('conflicting_duplicate_rejected',isempty(noAck) && numel(b.delivered)==1 && strcmp(event.name,'sequence_rejected'),'Same sequence with another CRC-valid payload is rejected');
[futureTx,~]=uaclink.transmit(uint8('future'),2,1,cfg);
[b,noAck,event]=uaclink.on_waveform(b,padded(futureTx,cfg),cfg,1.9);
rows{end+1}=test('unexpected_sequence_rejected',isempty(noAck) && numel(b.delivered)==1 && strcmp(event.name,'sequence_rejected'),'Expectedseq1 rejectsseq2');
[a,~,event]=uaclink.on_waveform(a,padded(ack,cfg),cfg,3.0);
rows{end+1}=test('matching_ack_clears_pending',strcmp(a.state,'idle') && a.next_sequence==1 && strcmp(event.name,'ack_received'),'CRC-valid matching acoustic ACK advances sequence');
% At the uint16 boundary both endpoints agree on their session counter.
a.next_sequence=65535; b.expected_sequence=65535;
[a,tx]=uaclink.send_message(a,'wrap',cfg,6);
[b,ack,~]=uaclink.on_waveform(b,padded(tx,cfg),cfg,7);
[a,~,~]=uaclink.on_waveform(a,padded(ack,cfg),cfg,8);
rows{end+1}=test('sequence_wrap',a.next_sequence==0 && b.expected_sequence==0 && strcmp(a.state,'idle'),'65535 wraps to0 for coordinated endpoints');
records=vertcat(rows{:}); writetable(struct2table(records),fullfile(out,'validation_cases.csv'));
report=struct('status','PASS','total_cases',numel(records),'passed_cases',sum([records.passed]), ...
    'all_passed',all([records.passed]),'crc_single_bit_trials',numel(bits), ...
    'noise_only_trials',10,'noise_only_false_accepts',falseAccept, ...
    'scope','Bounded software protocol + synthetic waveform tests only; no measured-channel or hardware claim', ...
    'receiver_inputs','received waveform and public cfg only; parsed length is decoded from64 header bits', ...
    'matlab_version',version,'created_utc',char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss''Z''')));
if ~report.all_passed, report.status='FAIL'; end
fid=fopen(fullfile(out,'validation.json'),'w'); assert(fid>=0); fwrite(fid,unicode2native(jsonencode(report),'UTF-8'),'uint8'); fclose(fid);
disp(struct2table(records));
fprintf('Packet protocol validation: %s (%d/%d)\n',report.status,report.passed_cases,report.total_cases);
assert(report.all_passed,'Packet protocol validation failed.');
end

function row=test(name,passed,details)
row=struct('name',name,'passed',logical(passed),'details',details);
end
function value=throws(action)
value=false; try, action(); catch, value=true; end
end
function rx=padded(tx,cfg)
rx=[zeros(round(.16*cfg.fs),1);tx;zeros(round(.12*cfg.fs),1)];
end
