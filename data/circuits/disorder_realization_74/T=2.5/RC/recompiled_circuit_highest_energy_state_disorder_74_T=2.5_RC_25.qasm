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
rz(-2.3251301) q[0];
sx q[0];
rz(-0.10179585) q[0];
sx q[0];
rz(2.6020004) q[0];
rz(-2.645283) q[1];
sx q[1];
rz(-2.8318383) q[1];
sx q[1];
rz(2.6024979) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2236299) q[0];
sx q[0];
rz(-1.5453891) q[0];
sx q[0];
rz(0.44141234) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7221872) q[2];
sx q[2];
rz(-1.9266204) q[2];
sx q[2];
rz(0.21499888) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2869563) q[1];
sx q[1];
rz(-2.7883734) q[1];
sx q[1];
rz(-1.9853206) q[1];
x q[2];
rz(-1.6786537) q[3];
sx q[3];
rz(-2.4898006) q[3];
sx q[3];
rz(0.45676132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3794136) q[2];
sx q[2];
rz(-2.2654686) q[2];
sx q[2];
rz(1.9525105) q[2];
rz(1.9573697) q[3];
sx q[3];
rz(-0.90457478) q[3];
sx q[3];
rz(-1.5949465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14520833) q[0];
sx q[0];
rz(-1.5518016) q[0];
sx q[0];
rz(1.3231963) q[0];
rz(-0.48201758) q[1];
sx q[1];
rz(-0.91164416) q[1];
sx q[1];
rz(0.97420305) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8708987) q[0];
sx q[0];
rz(-0.93241954) q[0];
sx q[0];
rz(-0.22298546) q[0];
rz(0.85136885) q[2];
sx q[2];
rz(-0.15002827) q[2];
sx q[2];
rz(-1.139037) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.43073359) q[1];
sx q[1];
rz(-1.3361592) q[1];
sx q[1];
rz(0.55107848) q[1];
x q[2];
rz(0.36688585) q[3];
sx q[3];
rz(-0.59795934) q[3];
sx q[3];
rz(0.46213996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.0049858967) q[2];
sx q[2];
rz(-0.58711457) q[2];
sx q[2];
rz(-0.74964398) q[2];
rz(-0.43198112) q[3];
sx q[3];
rz(-1.1155198) q[3];
sx q[3];
rz(-1.8384793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5169446) q[0];
sx q[0];
rz(-0.2722781) q[0];
sx q[0];
rz(0.6063478) q[0];
rz(-1.8079405) q[1];
sx q[1];
rz(-1.2140112) q[1];
sx q[1];
rz(-0.25951728) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5657046) q[0];
sx q[0];
rz(-1.4153) q[0];
sx q[0];
rz(-1.3208273) q[0];
rz(-pi) q[1];
rz(0.17510842) q[2];
sx q[2];
rz(-1.3710183) q[2];
sx q[2];
rz(1.4005043) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3644581) q[1];
sx q[1];
rz(-1.9455457) q[1];
sx q[1];
rz(2.0261637) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2058425) q[3];
sx q[3];
rz(-0.66422909) q[3];
sx q[3];
rz(0.85778824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2567265) q[2];
sx q[2];
rz(-1.7776411) q[2];
sx q[2];
rz(1.5941031) q[2];
rz(0.78440845) q[3];
sx q[3];
rz(-1.6268077) q[3];
sx q[3];
rz(-1.5756395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49482685) q[0];
sx q[0];
rz(-1.0972728) q[0];
sx q[0];
rz(2.8821017) q[0];
rz(1.9208113) q[1];
sx q[1];
rz(-2.6916598) q[1];
sx q[1];
rz(1.6993914) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4239765) q[0];
sx q[0];
rz(-1.0459131) q[0];
sx q[0];
rz(-2.5931326) q[0];
rz(-pi) q[1];
rz(0.13061009) q[2];
sx q[2];
rz(-2.161986) q[2];
sx q[2];
rz(-1.3863877) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4220548) q[1];
sx q[1];
rz(-1.9214848) q[1];
sx q[1];
rz(-1.8338051) q[1];
x q[2];
rz(-1.6569225) q[3];
sx q[3];
rz(-0.68400506) q[3];
sx q[3];
rz(0.67953426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1059025) q[2];
sx q[2];
rz(-1.8041958) q[2];
sx q[2];
rz(-0.40994677) q[2];
rz(-1.9299054) q[3];
sx q[3];
rz(-2.4544921) q[3];
sx q[3];
rz(1.80779) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.425151) q[0];
sx q[0];
rz(-0.71617675) q[0];
sx q[0];
rz(2.8675365) q[0];
rz(-0.43824276) q[1];
sx q[1];
rz(-1.2934877) q[1];
sx q[1];
rz(-2.4058707) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5490131) q[0];
sx q[0];
rz(-1.4837449) q[0];
sx q[0];
rz(2.3866231) q[0];
rz(-pi) q[1];
rz(-2.3724298) q[2];
sx q[2];
rz(-0.20649466) q[2];
sx q[2];
rz(2.1392876) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.085891) q[1];
sx q[1];
rz(-0.38905479) q[1];
sx q[1];
rz(-2.2772853) q[1];
rz(-2.2201331) q[3];
sx q[3];
rz(-1.9917352) q[3];
sx q[3];
rz(-2.4025847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.938544) q[2];
sx q[2];
rz(-1.6837348) q[2];
sx q[2];
rz(-0.61384821) q[2];
rz(2.7326873) q[3];
sx q[3];
rz(-2.1730065) q[3];
sx q[3];
rz(0.74234211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78855377) q[0];
sx q[0];
rz(-2.5034294) q[0];
sx q[0];
rz(-1.9045389) q[0];
rz(-1.4452112) q[1];
sx q[1];
rz(-0.45672363) q[1];
sx q[1];
rz(-0.17527418) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69578457) q[0];
sx q[0];
rz(-0.17477594) q[0];
sx q[0];
rz(-1.7150899) q[0];
rz(-pi) q[1];
rz(-1.0277599) q[2];
sx q[2];
rz(-2.1547085) q[2];
sx q[2];
rz(-1.2500545) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1239667) q[1];
sx q[1];
rz(-2.0793036) q[1];
sx q[1];
rz(1.0715436) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.82476576) q[3];
sx q[3];
rz(-2.3769393) q[3];
sx q[3];
rz(0.95461707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1767629) q[2];
sx q[2];
rz(-1.8076597) q[2];
sx q[2];
rz(2.1691587) q[2];
rz(-1.4771627) q[3];
sx q[3];
rz(-1.1494613) q[3];
sx q[3];
rz(-2.3798063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85457388) q[0];
sx q[0];
rz(-3.119097) q[0];
sx q[0];
rz(2.7884685) q[0];
rz(-2.1137721) q[1];
sx q[1];
rz(-1.9408344) q[1];
sx q[1];
rz(2.1305398) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7142657) q[0];
sx q[0];
rz(-1.2135226) q[0];
sx q[0];
rz(1.5861804) q[0];
rz(-0.89247668) q[2];
sx q[2];
rz(-0.8304285) q[2];
sx q[2];
rz(-0.76039808) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.87349975) q[1];
sx q[1];
rz(-2.756739) q[1];
sx q[1];
rz(1.1420239) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8500438) q[3];
sx q[3];
rz(-1.8729775) q[3];
sx q[3];
rz(2.9460689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8280243) q[2];
sx q[2];
rz(-0.32005969) q[2];
sx q[2];
rz(-1.3671406) q[2];
rz(1.8732871) q[3];
sx q[3];
rz(-1.6627848) q[3];
sx q[3];
rz(0.84793276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0276412) q[0];
sx q[0];
rz(-0.60751644) q[0];
sx q[0];
rz(-2.0116346) q[0];
rz(0.21895151) q[1];
sx q[1];
rz(-1.7349225) q[1];
sx q[1];
rz(-2.2311282) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6235891) q[0];
sx q[0];
rz(-1.9027038) q[0];
sx q[0];
rz(0.62923543) q[0];
x q[1];
rz(-1.6565645) q[2];
sx q[2];
rz(-1.0940486) q[2];
sx q[2];
rz(1.7976744) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7607019) q[1];
sx q[1];
rz(-2.5538951) q[1];
sx q[1];
rz(-2.0668849) q[1];
rz(-pi) q[2];
x q[2];
rz(0.45121737) q[3];
sx q[3];
rz(-2.3638551) q[3];
sx q[3];
rz(1.7672218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.8255446) q[2];
sx q[2];
rz(-2.3393708) q[2];
sx q[2];
rz(2.3160589) q[2];
rz(0.46142203) q[3];
sx q[3];
rz(-1.2666707) q[3];
sx q[3];
rz(-2.6958444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6398741) q[0];
sx q[0];
rz(-1.3383144) q[0];
sx q[0];
rz(-2.3097532) q[0];
rz(-0.72744751) q[1];
sx q[1];
rz(-1.7214382) q[1];
sx q[1];
rz(0.096253455) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.03503) q[0];
sx q[0];
rz(-1.0659395) q[0];
sx q[0];
rz(2.6120606) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5094069) q[2];
sx q[2];
rz(-1.824531) q[2];
sx q[2];
rz(0.77216567) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.61273328) q[1];
sx q[1];
rz(-1.7816356) q[1];
sx q[1];
rz(-0.092540578) q[1];
rz(-pi) q[2];
rz(2.2489266) q[3];
sx q[3];
rz(-0.24991194) q[3];
sx q[3];
rz(-1.753861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7353797) q[2];
sx q[2];
rz(-1.4344183) q[2];
sx q[2];
rz(1.5926788) q[2];
rz(2.5088572) q[3];
sx q[3];
rz(-1.2318719) q[3];
sx q[3];
rz(-2.8922141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(11/(9*pi)) q[0];
sx q[0];
rz(-1.0405552) q[0];
sx q[0];
rz(-2.7885875) q[0];
rz(0.32265916) q[1];
sx q[1];
rz(-0.32538515) q[1];
sx q[1];
rz(2.7696612) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3561365) q[0];
sx q[0];
rz(-1.8710941) q[0];
sx q[0];
rz(-2.9343283) q[0];
rz(-pi) q[1];
rz(-2.811889) q[2];
sx q[2];
rz(-1.7404544) q[2];
sx q[2];
rz(-0.28126954) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5185753) q[1];
sx q[1];
rz(-1.9037582) q[1];
sx q[1];
rz(-1.4013616) q[1];
x q[2];
rz(-1.6149893) q[3];
sx q[3];
rz(-0.72504504) q[3];
sx q[3];
rz(2.8738662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9670664) q[2];
sx q[2];
rz(-1.1897949) q[2];
sx q[2];
rz(-1.2971499) q[2];
rz(-0.8477115) q[3];
sx q[3];
rz(-1.1573557) q[3];
sx q[3];
rz(1.004809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90947718) q[0];
sx q[0];
rz(-1.3790601) q[0];
sx q[0];
rz(-2.7098304) q[0];
rz(-1.3879981) q[1];
sx q[1];
rz(-1.2139865) q[1];
sx q[1];
rz(-0.075275631) q[1];
rz(-2.3054988) q[2];
sx q[2];
rz(-2.2420364) q[2];
sx q[2];
rz(-1.4902761) q[2];
rz(2.344178) q[3];
sx q[3];
rz(-1.2684462) q[3];
sx q[3];
rz(1.0195337) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
