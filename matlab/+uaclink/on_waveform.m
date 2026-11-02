function [node, response, event] = on_waveform(node, waveform, cfg, nowS)
%ON_WAVEFORM Parse/CRC at each node; deliver once or emit actual modulated ACK.
response=zeros(0,1); packet=uaclink.receive(waveform,cfg);
event=struct('time_s',nowS,'name','rejected','packet_valid',packet.valid, ...
    'reason',packet.reason,'sequence',packet.sequence,'delivered',false);
if ~packet.valid, node.last_event='rejected'; return; end
if packet.type==2
    if strcmp(node.state,'waiting_ack') && packet.sequence==node.pending_sequence
        node.next_sequence=mod(node.pending_sequence+1,65536);
        node.pending_payload=zeros(1,0,'uint8'); node.pending_sequence=NaN;
        node.state='idle'; node.deadline_s=Inf; node.last_event='ack_received';
        event.name='ack_received'; event.reason='matching_ack';
    else
        node.last_event='unexpected_ack'; event.name='unexpected_ack';
        event.reason='no_matching_pending_sequence';
    end
    return;
end
if packet.sequence==node.expected_sequence
    node.delivered{end+1}=packet.payload;
    node.delivered_sequences(end+1)=packet.sequence;
    node.last_delivered_sequence=packet.sequence;
    node.last_delivered_payload=packet.payload;
    node.expected_sequence=mod(packet.sequence+1,65536);
    event.delivered=true; event.name='message_delivered'; event.reason='new_data';
elseif packet.sequence==node.last_delivered_sequence && ...
        isequal(packet.payload,node.last_delivered_payload)
    event.name='duplicate_acknowledged'; event.reason='duplicate_data';
else
    event.name='sequence_rejected'; event.reason='unexpected_or_conflicting_sequence';
    node.last_event=event.name; return;
end
[response,~]=uaclink.transmit(zeros(1,0,'uint8'),packet.sequence,2,cfg);
node.last_event=event.name;
end
