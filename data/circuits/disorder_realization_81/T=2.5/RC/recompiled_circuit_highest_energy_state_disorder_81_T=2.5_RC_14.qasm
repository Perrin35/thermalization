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
rz(-0.4001652) q[0];
sx q[0];
rz(-0.53181177) q[0];
sx q[0];
rz(2.3398633) q[0];
rz(0.62893686) q[1];
sx q[1];
rz(-1.8328272) q[1];
sx q[1];
rz(-2.7523249) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9297732) q[0];
sx q[0];
rz(-0.7302098) q[0];
sx q[0];
rz(-2.7279961) q[0];
rz(-pi) q[1];
x q[1];
rz(0.18125972) q[2];
sx q[2];
rz(-1.6158293) q[2];
sx q[2];
rz(-0.35895106) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.20922616) q[1];
sx q[1];
rz(-0.83360177) q[1];
sx q[1];
rz(0.7070138) q[1];
x q[2];
rz(1.5262768) q[3];
sx q[3];
rz(-2.4137437) q[3];
sx q[3];
rz(1.0893217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4929216) q[2];
sx q[2];
rz(-0.45009437) q[2];
sx q[2];
rz(2.8802803) q[2];
rz(2.785545) q[3];
sx q[3];
rz(-1.1084403) q[3];
sx q[3];
rz(-1.0935008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4556731) q[0];
sx q[0];
rz(-0.57722592) q[0];
sx q[0];
rz(-2.8976231) q[0];
rz(2.7120554) q[1];
sx q[1];
rz(-1.4725087) q[1];
sx q[1];
rz(2.4006749) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8888433) q[0];
sx q[0];
rz(-0.50562947) q[0];
sx q[0];
rz(-2.3028785) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6119611) q[2];
sx q[2];
rz(-1.1897693) q[2];
sx q[2];
rz(1.5059901) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5432905) q[1];
sx q[1];
rz(-1.0183698) q[1];
sx q[1];
rz(-1.8364947) q[1];
x q[2];
rz(-2.0617742) q[3];
sx q[3];
rz(-0.91697877) q[3];
sx q[3];
rz(-2.4308609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.7168768) q[2];
sx q[2];
rz(-1.5224785) q[2];
sx q[2];
rz(-0.87872046) q[2];
rz(2.1667513) q[3];
sx q[3];
rz(-1.2474371) q[3];
sx q[3];
rz(-1.6802906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9413302) q[0];
sx q[0];
rz(-0.55568475) q[0];
sx q[0];
rz(-0.19361198) q[0];
rz(-0.10784736) q[1];
sx q[1];
rz(-2.4272608) q[1];
sx q[1];
rz(0.81197369) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13511756) q[0];
sx q[0];
rz(-2.0251946) q[0];
sx q[0];
rz(-0.12250367) q[0];
rz(2.3247955) q[2];
sx q[2];
rz(-1.9765761) q[2];
sx q[2];
rz(-0.38015537) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1209315) q[1];
sx q[1];
rz(-1.963208) q[1];
sx q[1];
rz(-2.891178) q[1];
rz(1.1319086) q[3];
sx q[3];
rz(-1.4445856) q[3];
sx q[3];
rz(2.6708093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0783483) q[2];
sx q[2];
rz(-1.8599267) q[2];
sx q[2];
rz(0.0035302103) q[2];
rz(0.53623903) q[3];
sx q[3];
rz(-0.25501525) q[3];
sx q[3];
rz(0.53853881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2980767) q[0];
sx q[0];
rz(-1.902782) q[0];
sx q[0];
rz(-2.3115944) q[0];
rz(-1.4044546) q[1];
sx q[1];
rz(-2.2300215) q[1];
sx q[1];
rz(2.8703168) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1034086) q[0];
sx q[0];
rz(-1.8233946) q[0];
sx q[0];
rz(0.20471666) q[0];
x q[1];
rz(-3.0433367) q[2];
sx q[2];
rz(-2.4742804) q[2];
sx q[2];
rz(-1.7010393) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.30078706) q[1];
sx q[1];
rz(-1.5432165) q[1];
sx q[1];
rz(-0.3572398) q[1];
x q[2];
rz(1.4129407) q[3];
sx q[3];
rz(-0.78031555) q[3];
sx q[3];
rz(1.3387967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6259049) q[2];
sx q[2];
rz(-0.38334623) q[2];
sx q[2];
rz(2.3670727) q[2];
rz(-2.1382616) q[3];
sx q[3];
rz(-2.7480875) q[3];
sx q[3];
rz(-0.5459319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9970053) q[0];
sx q[0];
rz(-2.044675) q[0];
sx q[0];
rz(-1.0928094) q[0];
rz(-1.3358759) q[1];
sx q[1];
rz(-2.3922258) q[1];
sx q[1];
rz(-0.54378477) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8636759) q[0];
sx q[0];
rz(-1.5161235) q[0];
sx q[0];
rz(-1.3115505) q[0];
rz(-pi) q[1];
rz(-1.6751965) q[2];
sx q[2];
rz(-0.89492765) q[2];
sx q[2];
rz(0.75985786) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.136506) q[1];
sx q[1];
rz(-2.4509811) q[1];
sx q[1];
rz(2.4588799) q[1];
rz(-pi) q[2];
rz(0.001212515) q[3];
sx q[3];
rz(-2.1911484) q[3];
sx q[3];
rz(-1.9263489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.233923) q[2];
sx q[2];
rz(-2.449514) q[2];
sx q[2];
rz(-2.862759) q[2];
rz(-0.28576609) q[3];
sx q[3];
rz(-2.2279492) q[3];
sx q[3];
rz(0.41551503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7061507) q[0];
sx q[0];
rz(-0.64318648) q[0];
sx q[0];
rz(-2.1162794) q[0];
rz(-3.0030491) q[1];
sx q[1];
rz(-1.5899315) q[1];
sx q[1];
rz(2.9655546) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4885725) q[0];
sx q[0];
rz(-2.5551642) q[0];
sx q[0];
rz(-2.6676548) q[0];
rz(-pi) q[1];
x q[1];
rz(0.35961173) q[2];
sx q[2];
rz(-1.7533025) q[2];
sx q[2];
rz(2.3881641) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.306539) q[1];
sx q[1];
rz(-1.3471148) q[1];
sx q[1];
rz(2.064631) q[1];
rz(-1.995924) q[3];
sx q[3];
rz(-1.4151203) q[3];
sx q[3];
rz(-2.3666414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.90870086) q[2];
sx q[2];
rz(-0.42753926) q[2];
sx q[2];
rz(0.78988451) q[2];
rz(0.32957736) q[3];
sx q[3];
rz(-1.7695844) q[3];
sx q[3];
rz(-0.65143877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71094197) q[0];
sx q[0];
rz(-2.9281404) q[0];
sx q[0];
rz(-0.99367225) q[0];
rz(-1.0020533) q[1];
sx q[1];
rz(-2.1498945) q[1];
sx q[1];
rz(1.4659945) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77540175) q[0];
sx q[0];
rz(-0.68781536) q[0];
sx q[0];
rz(-2.8596876) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3411677) q[2];
sx q[2];
rz(-2.672019) q[2];
sx q[2];
rz(1.1498677) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8154303) q[1];
sx q[1];
rz(-0.23351352) q[1];
sx q[1];
rz(-0.68016078) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0577578) q[3];
sx q[3];
rz(-0.80082244) q[3];
sx q[3];
rz(1.458925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7953636) q[2];
sx q[2];
rz(-0.84285223) q[2];
sx q[2];
rz(0.25974926) q[2];
rz(-0.78884697) q[3];
sx q[3];
rz(-2.2987821) q[3];
sx q[3];
rz(1.747725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2271093) q[0];
sx q[0];
rz(-1.2863343) q[0];
sx q[0];
rz(2.1504543) q[0];
rz(1.8046851) q[1];
sx q[1];
rz(-2.3146345) q[1];
sx q[1];
rz(1.7327259) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0175655) q[0];
sx q[0];
rz(-0.47241898) q[0];
sx q[0];
rz(-0.25644619) q[0];
x q[1];
rz(-3.0524247) q[2];
sx q[2];
rz(-1.0341511) q[2];
sx q[2];
rz(-2.9111957) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0947021) q[1];
sx q[1];
rz(-0.45744236) q[1];
sx q[1];
rz(0.66384683) q[1];
rz(-pi) q[2];
rz(-1.4015607) q[3];
sx q[3];
rz(-2.2391161) q[3];
sx q[3];
rz(-3.075222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9764497) q[2];
sx q[2];
rz(-1.6565448) q[2];
sx q[2];
rz(2.5531947) q[2];
rz(-0.31383651) q[3];
sx q[3];
rz(-2.2321759) q[3];
sx q[3];
rz(-0.61217827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2345851) q[0];
sx q[0];
rz(-2.0186277) q[0];
sx q[0];
rz(-3.0871952) q[0];
rz(0.76039487) q[1];
sx q[1];
rz(-2.1957928) q[1];
sx q[1];
rz(2.0423582) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1791819) q[0];
sx q[0];
rz(-2.1405309) q[0];
sx q[0];
rz(0.197844) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.80087687) q[2];
sx q[2];
rz(-0.95716864) q[2];
sx q[2];
rz(-0.18979812) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1832378) q[1];
sx q[1];
rz(-2.3545165) q[1];
sx q[1];
rz(1.8182204) q[1];
x q[2];
rz(-0.051088574) q[3];
sx q[3];
rz(-1.6815485) q[3];
sx q[3];
rz(-2.1837019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2796563) q[2];
sx q[2];
rz(-1.6013689) q[2];
sx q[2];
rz(-2.3242365) q[2];
rz(3.127408) q[3];
sx q[3];
rz(-0.395702) q[3];
sx q[3];
rz(-2.0524249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31187439) q[0];
sx q[0];
rz(-2.0901966) q[0];
sx q[0];
rz(-2.9201087) q[0];
rz(2.7102176) q[1];
sx q[1];
rz(-0.90310496) q[1];
sx q[1];
rz(2.121714) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1843625) q[0];
sx q[0];
rz(-1.5910604) q[0];
sx q[0];
rz(3.094755) q[0];
x q[1];
rz(1.4657609) q[2];
sx q[2];
rz(-1.3702464) q[2];
sx q[2];
rz(2.5921164) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4012667) q[1];
sx q[1];
rz(-1.5366842) q[1];
sx q[1];
rz(1.9032065) q[1];
rz(-1.7027036) q[3];
sx q[3];
rz(-1.435809) q[3];
sx q[3];
rz(1.5832981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.43883103) q[2];
sx q[2];
rz(-1.2063382) q[2];
sx q[2];
rz(-1.8863401) q[2];
rz(0.025764763) q[3];
sx q[3];
rz(-1.8119101) q[3];
sx q[3];
rz(-2.6222031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41853607) q[0];
sx q[0];
rz(-2.3522455) q[0];
sx q[0];
rz(1.1242207) q[0];
rz(2.8070246) q[1];
sx q[1];
rz(-2.6523013) q[1];
sx q[1];
rz(2.1025067) q[1];
rz(-2.679297) q[2];
sx q[2];
rz(-1.8097327) q[2];
sx q[2];
rz(2.3859522) q[2];
rz(-2.1909621) q[3];
sx q[3];
rz(-1.9759569) q[3];
sx q[3];
rz(2.514545) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
