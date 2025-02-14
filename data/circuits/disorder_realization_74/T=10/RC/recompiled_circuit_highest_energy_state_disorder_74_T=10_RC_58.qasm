OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2008774) q[0];
sx q[0];
rz(-0.90056363) q[0];
sx q[0];
rz(-0.23298921) q[0];
rz(1.6147344) q[1];
sx q[1];
rz(-1.072071) q[1];
sx q[1];
rz(-1.1195247) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4716063) q[0];
sx q[0];
rz(-1.3363839) q[0];
sx q[0];
rz(-0.013201518) q[0];
rz(-0.018466516) q[2];
sx q[2];
rz(-2.5800843) q[2];
sx q[2];
rz(-0.36040053) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.2762766) q[1];
sx q[1];
rz(-2.7807693) q[1];
sx q[1];
rz(1.3555384) q[1];
rz(-0.4654765) q[3];
sx q[3];
rz(-0.66282907) q[3];
sx q[3];
rz(0.62191546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.69922525) q[2];
sx q[2];
rz(-2.0825601) q[2];
sx q[2];
rz(0.11288682) q[2];
rz(2.8804307) q[3];
sx q[3];
rz(-1.7966725) q[3];
sx q[3];
rz(-1.8908709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3498822) q[0];
sx q[0];
rz(-3.0630906) q[0];
sx q[0];
rz(3.0872524) q[0];
rz(0.21866523) q[1];
sx q[1];
rz(-1.668914) q[1];
sx q[1];
rz(-0.36453077) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4512149) q[0];
sx q[0];
rz(-1.0752819) q[0];
sx q[0];
rz(0.69083237) q[0];
rz(-pi) q[1];
x q[1];
rz(0.40016115) q[2];
sx q[2];
rz(-0.98099835) q[2];
sx q[2];
rz(1.9831744) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.042885429) q[1];
sx q[1];
rz(-1.646893) q[1];
sx q[1];
rz(-2.5185381) q[1];
x q[2];
rz(-1.909799) q[3];
sx q[3];
rz(-1.399013) q[3];
sx q[3];
rz(1.671227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.24078044) q[2];
sx q[2];
rz(-1.4419) q[2];
sx q[2];
rz(-1.1085294) q[2];
rz(2.945914) q[3];
sx q[3];
rz(-1.9340197) q[3];
sx q[3];
rz(1.6746563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7473258) q[0];
sx q[0];
rz(-2.1013923) q[0];
sx q[0];
rz(-2.105383) q[0];
rz(-0.58468435) q[1];
sx q[1];
rz(-1.6744637) q[1];
sx q[1];
rz(-2.5856957) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5591878) q[0];
sx q[0];
rz(-0.83600658) q[0];
sx q[0];
rz(2.56675) q[0];
x q[1];
rz(1.6481208) q[2];
sx q[2];
rz(-1.6082885) q[2];
sx q[2];
rz(2.897597) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.039274065) q[1];
sx q[1];
rz(-2.2688451) q[1];
sx q[1];
rz(-0.1711425) q[1];
rz(-pi) q[2];
rz(2.6494637) q[3];
sx q[3];
rz(-0.50072296) q[3];
sx q[3];
rz(0.89195097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3802152) q[2];
sx q[2];
rz(-2.9351202) q[2];
sx q[2];
rz(-1.9492487) q[2];
rz(-3.0377667) q[3];
sx q[3];
rz(-1.9793648) q[3];
sx q[3];
rz(1.9942358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1059859) q[0];
sx q[0];
rz(-2.5150531) q[0];
sx q[0];
rz(-1.4242127) q[0];
rz(0.72319889) q[1];
sx q[1];
rz(-0.23591787) q[1];
sx q[1];
rz(-0.22044388) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.047565) q[0];
sx q[0];
rz(-1.936543) q[0];
sx q[0];
rz(-2.0630552) q[0];
x q[1];
rz(2.1313558) q[2];
sx q[2];
rz(-1.3187485) q[2];
sx q[2];
rz(0.28365669) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4749149) q[1];
sx q[1];
rz(-1.4581469) q[1];
sx q[1];
rz(-1.9289609) q[1];
rz(0.8414874) q[3];
sx q[3];
rz(-0.53181767) q[3];
sx q[3];
rz(-1.8790661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7580938) q[2];
sx q[2];
rz(-1.8269962) q[2];
sx q[2];
rz(-0.51551762) q[2];
rz(-3.0564485) q[3];
sx q[3];
rz(-2.4862423) q[3];
sx q[3];
rz(-0.83318025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67007095) q[0];
sx q[0];
rz(-0.18297289) q[0];
sx q[0];
rz(0.66437379) q[0];
rz(-2.6145256) q[1];
sx q[1];
rz(-1.138569) q[1];
sx q[1];
rz(-1.1767496) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19163469) q[0];
sx q[0];
rz(-1.598888) q[0];
sx q[0];
rz(-1.4856443) q[0];
x q[1];
rz(-2.2559153) q[2];
sx q[2];
rz(-1.204245) q[2];
sx q[2];
rz(1.4057856) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.55498153) q[1];
sx q[1];
rz(-2.9025295) q[1];
sx q[1];
rz(-3.0596517) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.17950157) q[3];
sx q[3];
rz(-2.13509) q[3];
sx q[3];
rz(-1.5998942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1145757) q[2];
sx q[2];
rz(-2.9288374) q[2];
sx q[2];
rz(0.93811718) q[2];
rz(-1.045687) q[3];
sx q[3];
rz(-2.9595879) q[3];
sx q[3];
rz(-0.81030455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.806458) q[0];
sx q[0];
rz(-1.198575) q[0];
sx q[0];
rz(-1.2499811) q[0];
rz(0.20420034) q[1];
sx q[1];
rz(-0.67664346) q[1];
sx q[1];
rz(-2.2883889) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7366587) q[0];
sx q[0];
rz(-2.6879426) q[0];
sx q[0];
rz(-0.69866314) q[0];
rz(2.5218042) q[2];
sx q[2];
rz(-1.0519805) q[2];
sx q[2];
rz(-2.3960115) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0895929) q[1];
sx q[1];
rz(-1.9885364) q[1];
sx q[1];
rz(0.79302782) q[1];
rz(-pi) q[2];
rz(0.12746396) q[3];
sx q[3];
rz(-1.9098305) q[3];
sx q[3];
rz(1.6515428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0975254) q[2];
sx q[2];
rz(-1.7996457) q[2];
sx q[2];
rz(-1.1996783) q[2];
rz(-1.0058282) q[3];
sx q[3];
rz(-2.2483716) q[3];
sx q[3];
rz(-2.306166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47578874) q[0];
sx q[0];
rz(-2.8732193) q[0];
sx q[0];
rz(2.1589808) q[0];
rz(1.3081374) q[1];
sx q[1];
rz(-1.2160701) q[1];
sx q[1];
rz(1.7024202) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7076232) q[0];
sx q[0];
rz(-2.6585007) q[0];
sx q[0];
rz(2.3695495) q[0];
x q[1];
rz(-2.0928395) q[2];
sx q[2];
rz(-1.6380651) q[2];
sx q[2];
rz(-2.1362338) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5251617) q[1];
sx q[1];
rz(-2.3724864) q[1];
sx q[1];
rz(2.0517212) q[1];
x q[2];
rz(-1.4286707) q[3];
sx q[3];
rz(-2.0439337) q[3];
sx q[3];
rz(-2.9700235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.16284379) q[2];
sx q[2];
rz(-0.89357251) q[2];
sx q[2];
rz(-0.54235512) q[2];
rz(-2.1584623) q[3];
sx q[3];
rz(-1.1636846) q[3];
sx q[3];
rz(-2.2014309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5777609) q[0];
sx q[0];
rz(-1.7153808) q[0];
sx q[0];
rz(-2.3610709) q[0];
rz(-1.1753987) q[1];
sx q[1];
rz(-2.5927717) q[1];
sx q[1];
rz(-1.0343879) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7808409) q[0];
sx q[0];
rz(-0.409161) q[0];
sx q[0];
rz(2.7098814) q[0];
rz(-pi) q[1];
rz(-1.4277568) q[2];
sx q[2];
rz(-1.489822) q[2];
sx q[2];
rz(1.2317927) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3745034) q[1];
sx q[1];
rz(-1.4562291) q[1];
sx q[1];
rz(-2.7388045) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.54776056) q[3];
sx q[3];
rz(-1.8005074) q[3];
sx q[3];
rz(1.821777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.15022755) q[2];
sx q[2];
rz(-1.7346202) q[2];
sx q[2];
rz(-3.0312209) q[2];
rz(-1.2982093) q[3];
sx q[3];
rz(-1.7896264) q[3];
sx q[3];
rz(-2.8492294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89123911) q[0];
sx q[0];
rz(-2.1284916) q[0];
sx q[0];
rz(1.2497586) q[0];
rz(0.091014422) q[1];
sx q[1];
rz(-2.4200771) q[1];
sx q[1];
rz(-2.4299842) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34214333) q[0];
sx q[0];
rz(-0.79252386) q[0];
sx q[0];
rz(-1.5613129) q[0];
rz(-pi) q[1];
rz(0.72451024) q[2];
sx q[2];
rz(-0.9369999) q[2];
sx q[2];
rz(-2.8453635) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8636057) q[1];
sx q[1];
rz(-2.5558439) q[1];
sx q[1];
rz(-0.65237712) q[1];
rz(-pi) q[2];
x q[2];
rz(0.36104155) q[3];
sx q[3];
rz(-0.78892498) q[3];
sx q[3];
rz(1.5004683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.3332425) q[2];
sx q[2];
rz(-1.4834206) q[2];
sx q[2];
rz(1.1727772) q[2];
rz(-1.6138389) q[3];
sx q[3];
rz(-1.0781735) q[3];
sx q[3];
rz(-3.0465904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67977366) q[0];
sx q[0];
rz(-3.0175896) q[0];
sx q[0];
rz(-1.5661731) q[0];
rz(-2.7836986) q[1];
sx q[1];
rz(-2.0707668) q[1];
sx q[1];
rz(0.90248743) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47115147) q[0];
sx q[0];
rz(-2.4629399) q[0];
sx q[0];
rz(1.9891692) q[0];
x q[1];
rz(0.029377653) q[2];
sx q[2];
rz(-2.4973719) q[2];
sx q[2];
rz(2.6219896) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.984362) q[1];
sx q[1];
rz(-0.8547201) q[1];
sx q[1];
rz(-0.44992076) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4290513) q[3];
sx q[3];
rz(-0.55484164) q[3];
sx q[3];
rz(-2.6920163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1050528) q[2];
sx q[2];
rz(-0.5883216) q[2];
sx q[2];
rz(-1.6416637) q[2];
rz(-2.6203652) q[3];
sx q[3];
rz(-2.9977377) q[3];
sx q[3];
rz(2.1874793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70984107) q[0];
sx q[0];
rz(-2.3713645) q[0];
sx q[0];
rz(1.4159528) q[0];
rz(0.00090986666) q[1];
sx q[1];
rz(-1.4716499) q[1];
sx q[1];
rz(1.6843527) q[1];
rz(-0.098401423) q[2];
sx q[2];
rz(-0.99953166) q[2];
sx q[2];
rz(2.4649323) q[2];
rz(-2.3933181) q[3];
sx q[3];
rz(-1.0327485) q[3];
sx q[3];
rz(-2.2365981) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
