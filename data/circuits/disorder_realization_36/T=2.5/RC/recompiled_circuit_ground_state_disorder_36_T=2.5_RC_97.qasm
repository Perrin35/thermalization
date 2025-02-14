OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.26389709) q[0];
sx q[0];
rz(-1.0385624) q[0];
sx q[0];
rz(-0.39499083) q[0];
rz(0.78139961) q[1];
sx q[1];
rz(7.6548792) q[1];
sx q[1];
rz(10.21488) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.794644) q[0];
sx q[0];
rz(-2.145549) q[0];
sx q[0];
rz(-1.5331506) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0095135) q[2];
sx q[2];
rz(-0.84442511) q[2];
sx q[2];
rz(-1.4216636) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2826281) q[1];
sx q[1];
rz(-0.55759768) q[1];
sx q[1];
rz(-1.2931812) q[1];
rz(-1.155962) q[3];
sx q[3];
rz(-1.1006163) q[3];
sx q[3];
rz(3.1347203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.52749741) q[2];
sx q[2];
rz(-0.082753269) q[2];
sx q[2];
rz(2.5772074) q[2];
rz(2.2218521) q[3];
sx q[3];
rz(-2.4904833) q[3];
sx q[3];
rz(1.1434327) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98519242) q[0];
sx q[0];
rz(-2.1512478) q[0];
sx q[0];
rz(1.9957969) q[0];
rz(1.0076373) q[1];
sx q[1];
rz(-2.2538908) q[1];
sx q[1];
rz(-3.0767344) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50837171) q[0];
sx q[0];
rz(-2.457058) q[0];
sx q[0];
rz(0.9100911) q[0];
rz(-0.77936184) q[2];
sx q[2];
rz(-1.1032192) q[2];
sx q[2];
rz(3.0288966) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7520719) q[1];
sx q[1];
rz(-2.9557142) q[1];
sx q[1];
rz(3.0044772) q[1];
rz(-pi) q[2];
rz(0.81873881) q[3];
sx q[3];
rz(-2.4609339) q[3];
sx q[3];
rz(-0.37938947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2603944) q[2];
sx q[2];
rz(-3.03646) q[2];
sx q[2];
rz(2.6283188) q[2];
rz(-1.5567635) q[3];
sx q[3];
rz(-1.9624174) q[3];
sx q[3];
rz(-1.5796327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2960812) q[0];
sx q[0];
rz(-3.0873612) q[0];
sx q[0];
rz(0.84501141) q[0];
rz(0.71146479) q[1];
sx q[1];
rz(-0.80159801) q[1];
sx q[1];
rz(-0.62666384) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.959528) q[0];
sx q[0];
rz(-1.7668889) q[0];
sx q[0];
rz(-0.0053868731) q[0];
x q[1];
rz(-2.9128252) q[2];
sx q[2];
rz(-1.946515) q[2];
sx q[2];
rz(1.8919593) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.16594812) q[1];
sx q[1];
rz(-2.2774146) q[1];
sx q[1];
rz(-2.9099275) q[1];
rz(-pi) q[2];
rz(0.0073284433) q[3];
sx q[3];
rz(-1.0798608) q[3];
sx q[3];
rz(-2.1027264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6381548) q[2];
sx q[2];
rz(-1.2388836) q[2];
sx q[2];
rz(2.9729291) q[2];
rz(-0.1564129) q[3];
sx q[3];
rz(-0.63786879) q[3];
sx q[3];
rz(0.28353459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4446568) q[0];
sx q[0];
rz(-2.038027) q[0];
sx q[0];
rz(-2.7705833) q[0];
rz(-1.8095398) q[1];
sx q[1];
rz(-1.5947554) q[1];
sx q[1];
rz(0.36088774) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1754809) q[0];
sx q[0];
rz(-1.6663486) q[0];
sx q[0];
rz(0.086009228) q[0];
rz(-pi) q[1];
rz(-0.12207403) q[2];
sx q[2];
rz(-1.5549107) q[2];
sx q[2];
rz(-1.0779013) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3176705) q[1];
sx q[1];
rz(-0.63613632) q[1];
sx q[1];
rz(-2.0348861) q[1];
rz(-1.0293352) q[3];
sx q[3];
rz(-2.1730221) q[3];
sx q[3];
rz(0.15439872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0648301) q[2];
sx q[2];
rz(-0.37134376) q[2];
sx q[2];
rz(-1.4999464) q[2];
rz(1.0547981) q[3];
sx q[3];
rz(-1.0902371) q[3];
sx q[3];
rz(1.1928026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.763279) q[0];
sx q[0];
rz(-1.167647) q[0];
sx q[0];
rz(-1.8462697) q[0];
rz(3.0150705) q[1];
sx q[1];
rz(-2.449072) q[1];
sx q[1];
rz(1.893938) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5751094) q[0];
sx q[0];
rz(-1.2903443) q[0];
sx q[0];
rz(1.4202315) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3734002) q[2];
sx q[2];
rz(-1.8594701) q[2];
sx q[2];
rz(-1.4216258) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3646255) q[1];
sx q[1];
rz(-2.7086621) q[1];
sx q[1];
rz(0.51212279) q[1];
rz(-pi) q[2];
x q[2];
rz(0.99301569) q[3];
sx q[3];
rz(-0.3330001) q[3];
sx q[3];
rz(1.1525994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.36914545) q[2];
sx q[2];
rz(-0.36448604) q[2];
sx q[2];
rz(-0.44949731) q[2];
rz(-0.43826023) q[3];
sx q[3];
rz(-2.0354383) q[3];
sx q[3];
rz(-2.0800169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19675572) q[0];
sx q[0];
rz(-2.0773092) q[0];
sx q[0];
rz(-2.7979895) q[0];
rz(-0.65394872) q[1];
sx q[1];
rz(-2.3908354) q[1];
sx q[1];
rz(2.5508945) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9621679) q[0];
sx q[0];
rz(-2.5969167) q[0];
sx q[0];
rz(-1.1187798) q[0];
x q[1];
rz(-0.78738625) q[2];
sx q[2];
rz(-0.91286589) q[2];
sx q[2];
rz(-1.7869365) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.52378435) q[1];
sx q[1];
rz(-1.439938) q[1];
sx q[1];
rz(0.010556825) q[1];
x q[2];
rz(-0.88149397) q[3];
sx q[3];
rz(-2.1432264) q[3];
sx q[3];
rz(1.6400169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.58064738) q[2];
sx q[2];
rz(-2.2717768) q[2];
sx q[2];
rz(-0.84442863) q[2];
rz(2.1010418) q[3];
sx q[3];
rz(-1.2795762) q[3];
sx q[3];
rz(-1.9677264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3867253) q[0];
sx q[0];
rz(-0.52229053) q[0];
sx q[0];
rz(-2.0544384) q[0];
rz(0.48671752) q[1];
sx q[1];
rz(-1.646128) q[1];
sx q[1];
rz(-1.4901935) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33460245) q[0];
sx q[0];
rz(-0.86754543) q[0];
sx q[0];
rz(3.0038805) q[0];
rz(-pi) q[1];
rz(0.59972024) q[2];
sx q[2];
rz(-2.2968946) q[2];
sx q[2];
rz(-0.10657379) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9864038) q[1];
sx q[1];
rz(-1.1935368) q[1];
sx q[1];
rz(-2.6946486) q[1];
rz(0.20743117) q[3];
sx q[3];
rz(-1.279502) q[3];
sx q[3];
rz(-2.8136611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7955486) q[2];
sx q[2];
rz(-2.2717924) q[2];
sx q[2];
rz(0.33934936) q[2];
rz(0.36335534) q[3];
sx q[3];
rz(-2.8715869) q[3];
sx q[3];
rz(1.5615162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5245847) q[0];
sx q[0];
rz(-1.120485) q[0];
sx q[0];
rz(2.8560915) q[0];
rz(0.78954804) q[1];
sx q[1];
rz(-1.5317081) q[1];
sx q[1];
rz(-1.550536) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11131903) q[0];
sx q[0];
rz(-2.7919214) q[0];
sx q[0];
rz(-2.5969778) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8218846) q[2];
sx q[2];
rz(-2.6007915) q[2];
sx q[2];
rz(-0.021440949) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.68622491) q[1];
sx q[1];
rz(-1.1535201) q[1];
sx q[1];
rz(-2.5646006) q[1];
rz(-pi) q[2];
rz(-2.474349) q[3];
sx q[3];
rz(-1.5074456) q[3];
sx q[3];
rz(0.9982619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.32621128) q[2];
sx q[2];
rz(-2.9823494) q[2];
sx q[2];
rz(-0.86686575) q[2];
rz(0.73831144) q[3];
sx q[3];
rz(-1.5840931) q[3];
sx q[3];
rz(-2.2332634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0758122) q[0];
sx q[0];
rz(-2.3920238) q[0];
sx q[0];
rz(2.442389) q[0];
rz(0.54222822) q[1];
sx q[1];
rz(-1.7991853) q[1];
sx q[1];
rz(2.8453765) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64357054) q[0];
sx q[0];
rz(-1.6403461) q[0];
sx q[0];
rz(-2.3031917) q[0];
rz(-0.64546236) q[2];
sx q[2];
rz(-2.2156567) q[2];
sx q[2];
rz(-1.0114111) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.49668) q[1];
sx q[1];
rz(-1.0154004) q[1];
sx q[1];
rz(2.4128561) q[1];
x q[2];
rz(1.3557387) q[3];
sx q[3];
rz(-2.7557105) q[3];
sx q[3];
rz(-0.99487153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2518623) q[2];
sx q[2];
rz(-2.6984062) q[2];
sx q[2];
rz(2.6207793) q[2];
rz(0.03171799) q[3];
sx q[3];
rz(-0.40516502) q[3];
sx q[3];
rz(-3.090455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6601335) q[0];
sx q[0];
rz(-2.0429459) q[0];
sx q[0];
rz(-2.3302186) q[0];
rz(-0.42586455) q[1];
sx q[1];
rz(-0.90161294) q[1];
sx q[1];
rz(-1.3070377) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9872914) q[0];
sx q[0];
rz(-1.5736921) q[0];
sx q[0];
rz(-1.5737246) q[0];
rz(-pi) q[1];
rz(-3.1205253) q[2];
sx q[2];
rz(-2.0352007) q[2];
sx q[2];
rz(0.5201425) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9365434) q[1];
sx q[1];
rz(-1.2495511) q[1];
sx q[1];
rz(-0.46830362) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6622703) q[3];
sx q[3];
rz(-0.36277825) q[3];
sx q[3];
rz(-2.5517705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9797152) q[2];
sx q[2];
rz(-2.3780509) q[2];
sx q[2];
rz(-0.36583501) q[2];
rz(0.93519917) q[3];
sx q[3];
rz(-1.5949944) q[3];
sx q[3];
rz(0.38177761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3314421) q[0];
sx q[0];
rz(-1.7773542) q[0];
sx q[0];
rz(-1.6553028) q[0];
rz(-1.4797795) q[1];
sx q[1];
rz(-1.9098837) q[1];
sx q[1];
rz(2.2178537) q[1];
rz(-1.6681832) q[2];
sx q[2];
rz(-2.591572) q[2];
sx q[2];
rz(-0.77629065) q[2];
rz(-2.8769735) q[3];
sx q[3];
rz(-1.7195279) q[3];
sx q[3];
rz(-1.4715094) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
