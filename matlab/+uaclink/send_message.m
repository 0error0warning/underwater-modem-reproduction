function [node, waveform] = send_message(node, message, cfg, nowS)
%SEND_MESSAGE Queue one UTF-8 message; sender waits for a CRC-valid matching ACK.
assert(strcmp(node.state,'idle'), 'uaclink:SenderBusy', 'Sender is not idle.');
assert(ischar(message) || (isstring(message) && isscalar(message)), ...
    'uaclink:InvalidMessage', 'Expected one text message.');
payload=unicode2native(char(message),'UTF-8');
[waveform,info]=uaclink.transmit(payload,node.next_sequence,1,cfg);
node.pending_payload=uint8(payload); node.pending_sequence=node.next_sequence;
node.attempts=1; node.deadline_s=nowS+info.duration_s+cfg.link.ack_wait_s;
node.state='waiting_ack'; node.last_event='data_sent';
end
