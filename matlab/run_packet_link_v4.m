function report = run_packet_link_v4()
%RUN_PACKET_LINK_V4 Two software nodes exchanging DATA and acoustic-frame ACK.
% Deterministic synthetic test channels. Intentional loss is harness-only.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'matlab'));
out=fullfile(root,'results','v4','packet_link');
if ~exist(out,'dir'), mkdir(out); end
cfg=uaclink.default_config(); cfg.phase_rad=1.37;
message='你好，水声 Modem!';
scenarios={'normal','lose_first_data','lose_first_ack','lose_all_data','lose_all_ack'};
rows={}; events={};
for index=1:numel(scenarios)
    scenario=scenarios{index};
    a=uaclink.new_node('A'); b=uaclink.new_node('B'); nowS=0;
    [a,tx]=uaclink.send_message(a,message,cfg,nowS);
    while strcmp(a.state,'waiting_ack')
        attempt=a.attempts;
        event=struct('scenario',scenario,'attempt',attempt,'time_s',nowS, ...
            'data_reached_b',false,'data_crc_valid',false,'application_delivered',false, ...
            'ack_reached_a',false,'ack_crc_valid',false,'b_event','data_dropped', ...
            'a_event','waiting_ack','sequence',0);
        dropData=strcmp(scenario,'lose_all_data') || (strcmp(scenario,'lose_first_data') && attempt==1);
        if ~dropData
            rx=uaclink.synthetic_channel(tx,cfg,4000+index*100+attempt,0.5,20);
            arrivalS=nowS+numel(rx)/cfg.fs;
            [b,ack,bEvent]=uaclink.on_waveform(b,rx,cfg,arrivalS);
            event.data_reached_b=true; event.data_crc_valid=bEvent.packet_valid;
            event.application_delivered=bEvent.delivered; event.b_event=bEvent.name;
            dropAck=strcmp(scenario,'lose_all_ack') || (strcmp(scenario,'lose_first_ack') && attempt==1);
            if ~isempty(ack) && ~dropAck
                ackRx=uaclink.synthetic_channel(ack,cfg,8000+index*100+attempt,-0.5,20);
                ackArrivalS=arrivalS+numel(ackRx)/cfg.fs;
                assert(ackArrivalS<a.deadline_s,'ACK wait must exceed this test channel round trip.');
                [a,~,aEvent]=uaclink.on_waveform(a,ackRx,cfg,ackArrivalS);
                event.ack_reached_a=true; event.ack_crc_valid=aEvent.packet_valid;
                event.a_event=aEvent.name;
            end
        end
        events{end+1}=event; %#ok<AGROW>
        if strcmp(a.state,'waiting_ack')
            nowS=a.deadline_s;
            [a,tx]=uaclink.on_timeout(a,cfg,nowS);
        end
    end
    deliveredText='';
    if ~isempty(b.delivered), deliveredText=native2unicode(b.delivered{1},'UTF-8'); end
    expectedAttempts=1;
    if ismember(scenario,{'lose_first_data','lose_first_ack'}), expectedAttempts=2; end
    if ismember(scenario,{'lose_all_data','lose_all_ack'}), expectedAttempts=cfg.link.max_attempts; end
    expectedDelivery=~strcmp(scenario,'lose_all_data');
    expectedConfirmed=~ismember(scenario,{'lose_all_data','lose_all_ack'});
    passed=a.attempts==expectedAttempts && numel(b.delivered)==expectedDelivery && ...
        strcmp(a.state,'idle')==expectedConfirmed && ...
        (~expectedDelivery || strcmp(message,deliveredText));
    rows{end+1}=struct('scenario',scenario,'attempts',a.attempts, ...
        'sender_state',a.state,'sender_confirmed',strcmp(a.state,'idle'), ...
        'receiver_deliveries',numel(b.delivered),'received_text',deliveredText,'message_equal',strcmp(message,deliveredText), ...
        'passed',passed); %#ok<AGROW>
    fprintf('%s: attempts=%d state=%s deliveries=%d pass=%d\n',scenario,a.attempts,a.state,numel(b.delivered),passed);
    if ~isempty(deliveredText), fprintf('Node B output: %s\n',deliveredText); end
end
writetable(struct2table(vertcat(rows{:})),fullfile(out,'demo_scenarios.csv'));
writetable(struct2table(vertcat(events{:})),fullfile(out,'demo_events.csv'));
[exampleWire,~]=uaclink.encode_frame(unicode2native(message,'UTF-8'),0,1,cfg);
report=struct('status','PASS','scope','Synthetic delayed/Doppler/AWGN waveform DATA and ACK; not measured-channel or hardware validation', ...
    'message',message,'config',cfg,'scenario_count',numel(rows),'event_count',numel(events), ...
    'all_passed',all(cellfun(@(r)r.passed,rows)),'scenarios',{rows},'example_wire_bytes',exampleWire, ...
    'truth_isolation','Node APIs see only rx waveform, public cfg and local time; loss and physical parameters are harness-only.', ...
    'exhaustion_semantics','Delivery remains unconfirmed at A; B may have already delivered when all ACKs are lost.', ...
    'matlab_version',version,'created_utc',char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss''Z''')));
if ~report.all_passed, report.status='FAIL'; end
fid=fopen(fullfile(out,'demo.json'),'w'); assert(fid>=0); fwrite(fid,unicode2native(jsonencode(report),'UTF-8'),'uint8'); fclose(fid);
assert(report.all_passed,'Packet link demo failed.');
end
