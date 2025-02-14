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
rz(4.0539157) q[0];
sx q[0];
rz(11.231448) q[0];
rz(-1.1420684) q[1];
sx q[1];
rz(-0.51841441) q[1];
sx q[1];
rz(-1.4919182) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8700884) q[0];
sx q[0];
rz(-0.056500204) q[0];
sx q[0];
rz(1.2369878) q[0];
rz(0.71662997) q[2];
sx q[2];
rz(-2.306005) q[2];
sx q[2];
rz(0.40975964) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8644997) q[1];
sx q[1];
rz(-0.96597176) q[1];
sx q[1];
rz(-2.7983309) q[1];
rz(1.8133477) q[3];
sx q[3];
rz(-1.6036766) q[3];
sx q[3];
rz(1.3803079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7923183) q[2];
sx q[2];
rz(-0.26733843) q[2];
sx q[2];
rz(-3.086669) q[2];
rz(1.6867636) q[3];
sx q[3];
rz(-2.7456386) q[3];
sx q[3];
rz(1.5168064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9216264) q[0];
sx q[0];
rz(-1.3082137) q[0];
sx q[0];
rz(0.24800214) q[0];
rz(1.2451046) q[1];
sx q[1];
rz(-2.4047132) q[1];
sx q[1];
rz(-1.2287593) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1708536) q[0];
sx q[0];
rz(-1.5376933) q[0];
sx q[0];
rz(-2.956361) q[0];
rz(-pi) q[1];
rz(2.3360148) q[2];
sx q[2];
rz(-2.4586939) q[2];
sx q[2];
rz(2.3685092) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.23763021) q[1];
sx q[1];
rz(-1.7824689) q[1];
sx q[1];
rz(-0.71118506) q[1];
rz(-pi) q[2];
rz(1.5108438) q[3];
sx q[3];
rz(-1.9904117) q[3];
sx q[3];
rz(2.7762716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2174786) q[2];
sx q[2];
rz(-0.83196297) q[2];
sx q[2];
rz(2.5577616) q[2];
rz(-2.7122688) q[3];
sx q[3];
rz(-1.2238294) q[3];
sx q[3];
rz(2.1302285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9871224) q[0];
sx q[0];
rz(-0.38726375) q[0];
sx q[0];
rz(1.467147) q[0];
rz(1.6636498) q[1];
sx q[1];
rz(-1.4091622) q[1];
sx q[1];
rz(-2.0416226) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8505658) q[0];
sx q[0];
rz(-1.675559) q[0];
sx q[0];
rz(0.19922231) q[0];
rz(0.28950615) q[2];
sx q[2];
rz(-1.2443674) q[2];
sx q[2];
rz(3.141186) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5051859) q[1];
sx q[1];
rz(-2.6224394) q[1];
sx q[1];
rz(2.4207508) q[1];
x q[2];
rz(-0.44902841) q[3];
sx q[3];
rz(-2.0981023) q[3];
sx q[3];
rz(2.9849986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.786342) q[2];
sx q[2];
rz(-2.2189271) q[2];
sx q[2];
rz(1.5477017) q[2];
rz(1.330438) q[3];
sx q[3];
rz(-1.3733613) q[3];
sx q[3];
rz(1.0734585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7610157) q[0];
sx q[0];
rz(-1.3207734) q[0];
sx q[0];
rz(2.9837578) q[0];
rz(-0.24179587) q[1];
sx q[1];
rz(-0.38547412) q[1];
sx q[1];
rz(2.6604624) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0711827) q[0];
sx q[0];
rz(-2.2803223) q[0];
sx q[0];
rz(1.6874403) q[0];
rz(-pi) q[1];
rz(2.535216) q[2];
sx q[2];
rz(-2.7131792) q[2];
sx q[2];
rz(-2.1386752) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.85192902) q[1];
sx q[1];
rz(-0.73411513) q[1];
sx q[1];
rz(2.0925702) q[1];
rz(-2.1948482) q[3];
sx q[3];
rz(-0.3487772) q[3];
sx q[3];
rz(2.6117976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2699997) q[2];
sx q[2];
rz(-2.3303878) q[2];
sx q[2];
rz(0.2743741) q[2];
rz(1.6756049) q[3];
sx q[3];
rz(-1.2803187) q[3];
sx q[3];
rz(-2.2082641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6896553) q[0];
sx q[0];
rz(-0.87757293) q[0];
sx q[0];
rz(-2.080132) q[0];
rz(-0.44867107) q[1];
sx q[1];
rz(-0.45495382) q[1];
sx q[1];
rz(1.4400858) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7730624) q[0];
sx q[0];
rz(-2.321918) q[0];
sx q[0];
rz(1.3557662) q[0];
rz(-pi) q[1];
rz(-3.0147047) q[2];
sx q[2];
rz(-1.1694307) q[2];
sx q[2];
rz(-3.1201862) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.49855194) q[1];
sx q[1];
rz(-1.3128377) q[1];
sx q[1];
rz(-0.027679701) q[1];
rz(2.1813356) q[3];
sx q[3];
rz(-1.4818076) q[3];
sx q[3];
rz(-1.2582092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7777286) q[2];
sx q[2];
rz(-1.2960459) q[2];
sx q[2];
rz(1.9484005) q[2];
rz(-2.7423972) q[3];
sx q[3];
rz(-1.4123071) q[3];
sx q[3];
rz(-3.1138368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6804009) q[0];
sx q[0];
rz(-1.8541279) q[0];
sx q[0];
rz(-2.4611018) q[0];
rz(0.83241278) q[1];
sx q[1];
rz(-1.8871555) q[1];
sx q[1];
rz(-0.15636538) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.098861) q[0];
sx q[0];
rz(-1.8462688) q[0];
sx q[0];
rz(-1.5163438) q[0];
rz(0.83145468) q[2];
sx q[2];
rz(-1.0359087) q[2];
sx q[2];
rz(1.8597459) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.27626565) q[1];
sx q[1];
rz(-0.4368096) q[1];
sx q[1];
rz(0.86658367) q[1];
rz(-pi) q[2];
x q[2];
rz(0.12229192) q[3];
sx q[3];
rz(-1.7867136) q[3];
sx q[3];
rz(-1.9401039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.51123315) q[2];
sx q[2];
rz(-0.66143051) q[2];
sx q[2];
rz(0.88999256) q[2];
rz(0.47641274) q[3];
sx q[3];
rz(-2.2178631) q[3];
sx q[3];
rz(0.43058968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7156859) q[0];
sx q[0];
rz(-1.8222734) q[0];
sx q[0];
rz(2.529378) q[0];
rz(-1.3587492) q[1];
sx q[1];
rz(-0.76018676) q[1];
sx q[1];
rz(-2.8403958) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33443794) q[0];
sx q[0];
rz(-2.8557192) q[0];
sx q[0];
rz(1.8854333) q[0];
rz(-2.7686504) q[2];
sx q[2];
rz(-2.1939447) q[2];
sx q[2];
rz(-2.9187893) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0665022) q[1];
sx q[1];
rz(-1.1937703) q[1];
sx q[1];
rz(-2.6860363) q[1];
rz(-pi) q[2];
rz(1.8101569) q[3];
sx q[3];
rz(-1.3151957) q[3];
sx q[3];
rz(1.1472697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.907454) q[2];
sx q[2];
rz(-0.80628866) q[2];
sx q[2];
rz(-1.0478728) q[2];
rz(0.35510865) q[3];
sx q[3];
rz(-2.5706048) q[3];
sx q[3];
rz(0.6855489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2316786) q[0];
sx q[0];
rz(-2.7993918) q[0];
sx q[0];
rz(2.8177596) q[0];
rz(0.58770761) q[1];
sx q[1];
rz(-2.4617742) q[1];
sx q[1];
rz(-0.30276611) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.628474) q[0];
sx q[0];
rz(-1.895005) q[0];
sx q[0];
rz(-1.3036435) q[0];
rz(-pi) q[1];
rz(0.19925929) q[2];
sx q[2];
rz(-2.8997921) q[2];
sx q[2];
rz(-1.3346145) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7851315) q[1];
sx q[1];
rz(-1.0279127) q[1];
sx q[1];
rz(-0.42713844) q[1];
rz(1.15074) q[3];
sx q[3];
rz(-1.8057293) q[3];
sx q[3];
rz(0.7701503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6451463) q[2];
sx q[2];
rz(-2.1108997) q[2];
sx q[2];
rz(-2.5808064) q[2];
rz(2.136611) q[3];
sx q[3];
rz(-1.2882261) q[3];
sx q[3];
rz(2.0463478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0305369) q[0];
sx q[0];
rz(-0.66937864) q[0];
sx q[0];
rz(-2.3916767) q[0];
rz(0.38052446) q[1];
sx q[1];
rz(-2.1546202) q[1];
sx q[1];
rz(2.1662625) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0935287) q[0];
sx q[0];
rz(-0.87665999) q[0];
sx q[0];
rz(0.82406901) q[0];
x q[1];
rz(1.2911694) q[2];
sx q[2];
rz(-1.6284918) q[2];
sx q[2];
rz(-1.9967494) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1854461) q[1];
sx q[1];
rz(-2.5393894) q[1];
sx q[1];
rz(2.5700997) q[1];
rz(-pi) q[2];
rz(2.1998243) q[3];
sx q[3];
rz(-0.36702752) q[3];
sx q[3];
rz(1.1197156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.19899496) q[2];
sx q[2];
rz(-2.2364605) q[2];
sx q[2];
rz(-2.060176) q[2];
rz(3.0723451) q[3];
sx q[3];
rz(-1.705575) q[3];
sx q[3];
rz(-2.2190905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0104495) q[0];
sx q[0];
rz(-0.32587019) q[0];
sx q[0];
rz(0.3057873) q[0];
rz(-0.30793134) q[1];
sx q[1];
rz(-1.4762069) q[1];
sx q[1];
rz(-1.1526795) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9891805) q[0];
sx q[0];
rz(-2.3539641) q[0];
sx q[0];
rz(0.37352011) q[0];
x q[1];
rz(-0.20736097) q[2];
sx q[2];
rz(-2.2054858) q[2];
sx q[2];
rz(-2.7864252) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4977485) q[1];
sx q[1];
rz(-1.9744999) q[1];
sx q[1];
rz(-0.35404842) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3373371) q[3];
sx q[3];
rz(-0.065107927) q[3];
sx q[3];
rz(2.5447735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.61047381) q[2];
sx q[2];
rz(-1.0627154) q[2];
sx q[2];
rz(-1.4531685) q[2];
rz(0.48062634) q[3];
sx q[3];
rz(-2.5745945) q[3];
sx q[3];
rz(1.7467197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4881445) q[0];
sx q[0];
rz(-1.5339889) q[0];
sx q[0];
rz(1.4657159) q[0];
rz(2.0987971) q[1];
sx q[1];
rz(-1.7386309) q[1];
sx q[1];
rz(2.244619) q[1];
rz(-1.8389134) q[2];
sx q[2];
rz(-1.0514048) q[2];
sx q[2];
rz(1.9197293) q[2];
rz(-3.0328209) q[3];
sx q[3];
rz(-0.36009195) q[3];
sx q[3];
rz(2.910955) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
