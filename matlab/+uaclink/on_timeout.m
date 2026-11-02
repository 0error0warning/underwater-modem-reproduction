function [node, waveform] = on_timeout(node, cfg, nowS)
%ON_TIMEOUT Retransmit after deadline; exhausted means delivery is UNCONFIRMED.
waveform=zeros(0,1);
if ~strcmp(node.state,'waiting_ack') || nowS < node.deadline_s, return; end
if node.attempts >= cfg.link.max_attempts
    node.state='exhausted'; node.last_event='delivery_unconfirmed';
    node.deadline_s=Inf; return;
end
[waveform,info]=uaclink.transmit(node.pending_payload,node.pending_sequence,1,cfg);
node.attempts=node.attempts+1;
node.deadline_s=nowS+info.duration_s+cfg.link.ack_wait_s;
node.last_event='data_retransmitted';
end
