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
rz(0.063989446) q[0];
sx q[0];
rz(-2.2292697) q[0];
sx q[0];
rz(1.8066701) q[0];
rz(-1.1420684) q[1];
sx q[1];
rz(-0.51841441) q[1];
sx q[1];
rz(-1.4919182) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27150422) q[0];
sx q[0];
rz(-0.056500204) q[0];
sx q[0];
rz(1.2369878) q[0];
rz(-pi) q[1];
rz(-0.71662997) q[2];
sx q[2];
rz(-0.83558768) q[2];
sx q[2];
rz(0.40975964) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.277093) q[1];
sx q[1];
rz(-2.1756209) q[1];
sx q[1];
rz(2.7983309) q[1];
rz(-pi) q[2];
x q[2];
rz(1.434695) q[3];
sx q[3];
rz(-0.24472642) q[3];
sx q[3];
rz(2.8190106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7923183) q[2];
sx q[2];
rz(-2.8742542) q[2];
sx q[2];
rz(0.054923687) q[2];
rz(-1.4548291) q[3];
sx q[3];
rz(-0.39595404) q[3];
sx q[3];
rz(1.6247862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9216264) q[0];
sx q[0];
rz(-1.3082137) q[0];
sx q[0];
rz(-2.8935905) q[0];
rz(-1.8964881) q[1];
sx q[1];
rz(-2.4047132) q[1];
sx q[1];
rz(1.9128333) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7354483) q[0];
sx q[0];
rz(-1.3856674) q[0];
sx q[0];
rz(-1.6044751) q[0];
x q[1];
rz(2.6284559) q[2];
sx q[2];
rz(-2.0433132) q[2];
sx q[2];
rz(-3.0233011) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5723443) q[1];
sx q[1];
rz(-2.4048785) q[1];
sx q[1];
rz(0.31804292) q[1];
x q[2];
rz(-3.008083) q[3];
sx q[3];
rz(-2.7179681) q[3];
sx q[3];
rz(2.6299919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9241141) q[2];
sx q[2];
rz(-2.3096297) q[2];
sx q[2];
rz(2.5577616) q[2];
rz(0.42932388) q[3];
sx q[3];
rz(-1.9177633) q[3];
sx q[3];
rz(1.0113641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15447021) q[0];
sx q[0];
rz(-2.7543289) q[0];
sx q[0];
rz(1.6744457) q[0];
rz(1.6636498) q[1];
sx q[1];
rz(-1.7324305) q[1];
sx q[1];
rz(-1.0999701) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29102688) q[0];
sx q[0];
rz(-1.675559) q[0];
sx q[0];
rz(0.19922231) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9103545) q[2];
sx q[2];
rz(-1.2969839) q[2];
sx q[2];
rz(1.4751612) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5051859) q[1];
sx q[1];
rz(-0.51915324) q[1];
sx q[1];
rz(-2.4207508) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1446225) q[3];
sx q[3];
rz(-1.9553292) q[3];
sx q[3];
rz(-1.1763193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.786342) q[2];
sx q[2];
rz(-0.92266551) q[2];
sx q[2];
rz(-1.5938909) q[2];
rz(1.8111546) q[3];
sx q[3];
rz(-1.7682313) q[3];
sx q[3];
rz(-2.0681341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.380577) q[0];
sx q[0];
rz(-1.3207734) q[0];
sx q[0];
rz(2.9837578) q[0];
rz(2.8997968) q[1];
sx q[1];
rz(-0.38547412) q[1];
sx q[1];
rz(2.6604624) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24836981) q[0];
sx q[0];
rz(-2.4241894) q[0];
sx q[0];
rz(-0.1347085) q[0];
rz(-pi) q[1];
x q[1];
rz(2.535216) q[2];
sx q[2];
rz(-2.7131792) q[2];
sx q[2];
rz(-2.1386752) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.85192902) q[1];
sx q[1];
rz(-2.4074775) q[1];
sx q[1];
rz(-1.0490225) q[1];
rz(1.8577544) q[3];
sx q[3];
rz(-1.7718414) q[3];
sx q[3];
rz(-1.505681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2699997) q[2];
sx q[2];
rz(-2.3303878) q[2];
sx q[2];
rz(0.2743741) q[2];
rz(-1.4659878) q[3];
sx q[3];
rz(-1.2803187) q[3];
sx q[3];
rz(-2.2082641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45193732) q[0];
sx q[0];
rz(-2.2640197) q[0];
sx q[0];
rz(2.080132) q[0];
rz(-2.6929216) q[1];
sx q[1];
rz(-0.45495382) q[1];
sx q[1];
rz(-1.4400858) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0827918) q[0];
sx q[0];
rz(-0.77547204) q[0];
sx q[0];
rz(2.916921) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8605804) q[2];
sx q[2];
rz(-2.7216879) q[2];
sx q[2];
rz(0.29422255) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.49855194) q[1];
sx q[1];
rz(-1.3128377) q[1];
sx q[1];
rz(-3.113913) q[1];
rz(-2.1813356) q[3];
sx q[3];
rz(-1.4818076) q[3];
sx q[3];
rz(-1.8833835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3638641) q[2];
sx q[2];
rz(-1.2960459) q[2];
sx q[2];
rz(1.9484005) q[2];
rz(-0.39919546) q[3];
sx q[3];
rz(-1.7292855) q[3];
sx q[3];
rz(0.027755888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4611918) q[0];
sx q[0];
rz(-1.2874648) q[0];
sx q[0];
rz(2.4611018) q[0];
rz(-2.3091799) q[1];
sx q[1];
rz(-1.8871555) q[1];
sx q[1];
rz(2.9852273) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.844961) q[0];
sx q[0];
rz(-0.28066844) q[0];
sx q[0];
rz(-0.190221) q[0];
rz(-pi) q[1];
x q[1];
rz(0.84951289) q[2];
sx q[2];
rz(-2.2597183) q[2];
sx q[2];
rz(-2.920814) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1121516) q[1];
sx q[1];
rz(-1.2425155) q[1];
sx q[1];
rz(0.29354696) q[1];
rz(-pi) q[2];
rz(-1.353305) q[3];
sx q[3];
rz(-1.6902349) q[3];
sx q[3];
rz(-2.7986106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6303595) q[2];
sx q[2];
rz(-2.4801621) q[2];
sx q[2];
rz(-2.2516001) q[2];
rz(2.6651799) q[3];
sx q[3];
rz(-0.92372957) q[3];
sx q[3];
rz(0.43058968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42590672) q[0];
sx q[0];
rz(-1.3193193) q[0];
sx q[0];
rz(0.61221468) q[0];
rz(1.3587492) q[1];
sx q[1];
rz(-0.76018676) q[1];
sx q[1];
rz(-0.30119687) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66147236) q[0];
sx q[0];
rz(-1.8422707) q[0];
sx q[0];
rz(-0.090711509) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7686504) q[2];
sx q[2];
rz(-2.1939447) q[2];
sx q[2];
rz(-2.9187893) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4674356) q[1];
sx q[1];
rz(-1.1493719) q[1];
sx q[1];
rz(-1.9860877) q[1];
rz(-pi) q[2];
x q[2];
rz(2.878827) q[3];
sx q[3];
rz(-1.802236) q[3];
sx q[3];
rz(0.48514807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.23413868) q[2];
sx q[2];
rz(-2.335304) q[2];
sx q[2];
rz(-2.0937199) q[2];
rz(-0.35510865) q[3];
sx q[3];
rz(-0.57098782) q[3];
sx q[3];
rz(-2.4560438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.909914) q[0];
sx q[0];
rz(-2.7993918) q[0];
sx q[0];
rz(-0.32383305) q[0];
rz(0.58770761) q[1];
sx q[1];
rz(-2.4617742) q[1];
sx q[1];
rz(2.8388265) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2229109) q[0];
sx q[0];
rz(-0.41707539) q[0];
sx q[0];
rz(-2.4757451) q[0];
rz(1.6195756) q[2];
sx q[2];
rz(-1.80772) q[2];
sx q[2];
rz(1.1295484) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.3564612) q[1];
sx q[1];
rz(-2.11368) q[1];
sx q[1];
rz(-0.42713844) q[1];
x q[2];
rz(-1.9908526) q[3];
sx q[3];
rz(-1.3358634) q[3];
sx q[3];
rz(-0.7701503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4964464) q[2];
sx q[2];
rz(-1.030693) q[2];
sx q[2];
rz(0.56078625) q[2];
rz(-2.136611) q[3];
sx q[3];
rz(-1.2882261) q[3];
sx q[3];
rz(1.0952449) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1110558) q[0];
sx q[0];
rz(-0.66937864) q[0];
sx q[0];
rz(0.74991599) q[0];
rz(-2.7610682) q[1];
sx q[1];
rz(-0.98697248) q[1];
sx q[1];
rz(0.97533018) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98816865) q[0];
sx q[0];
rz(-1.0215217) q[0];
sx q[0];
rz(-0.84814056) q[0];
rz(-pi) q[1];
rz(-3.0815711) q[2];
sx q[2];
rz(-1.8499455) q[2];
sx q[2];
rz(2.699083) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.9561466) q[1];
sx q[1];
rz(-2.5393894) q[1];
sx q[1];
rz(2.5700997) q[1];
rz(-2.1998243) q[3];
sx q[3];
rz(-2.7745651) q[3];
sx q[3];
rz(1.1197156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9425977) q[2];
sx q[2];
rz(-0.90513217) q[2];
sx q[2];
rz(-1.0814166) q[2];
rz(-0.069247581) q[3];
sx q[3];
rz(-1.705575) q[3];
sx q[3];
rz(-2.2190905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0104495) q[0];
sx q[0];
rz(-0.32587019) q[0];
sx q[0];
rz(-2.8358054) q[0];
rz(0.30793134) q[1];
sx q[1];
rz(-1.4762069) q[1];
sx q[1];
rz(-1.9889132) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.688153) q[0];
sx q[0];
rz(-1.3092293) q[0];
sx q[0];
rz(0.75193172) q[0];
x q[1];
rz(-0.20736097) q[2];
sx q[2];
rz(-0.93610686) q[2];
sx q[2];
rz(2.7864252) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4977485) q[1];
sx q[1];
rz(-1.1670928) q[1];
sx q[1];
rz(0.35404842) q[1];
rz(-pi) q[2];
rz(1.5074499) q[3];
sx q[3];
rz(-1.5858486) q[3];
sx q[3];
rz(-1.2069595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5311188) q[2];
sx q[2];
rz(-2.0788772) q[2];
sx q[2];
rz(1.6884241) q[2];
rz(-0.48062634) q[3];
sx q[3];
rz(-0.56699816) q[3];
sx q[3];
rz(1.7467197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4881445) q[0];
sx q[0];
rz(-1.5339889) q[0];
sx q[0];
rz(1.4657159) q[0];
rz(1.0427955) q[1];
sx q[1];
rz(-1.4029618) q[1];
sx q[1];
rz(-0.89697368) q[1];
rz(-2.6063812) q[2];
sx q[2];
rz(-1.3387398) q[2];
sx q[2];
rz(-2.6571318) q[2];
rz(0.10877175) q[3];
sx q[3];
rz(-0.36009195) q[3];
sx q[3];
rz(2.910955) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
