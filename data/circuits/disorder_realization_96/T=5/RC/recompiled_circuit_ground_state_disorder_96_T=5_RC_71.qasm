OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.191303) q[0];
sx q[0];
rz(-2.8714955) q[0];
sx q[0];
rz(-0.88859963) q[0];
rz(1.8285881) q[1];
sx q[1];
rz(4.7409952) q[1];
sx q[1];
rz(10.805605) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61533538) q[0];
sx q[0];
rz(-1.3852296) q[0];
sx q[0];
rz(2.944817) q[0];
rz(-pi) q[1];
rz(-2.6216952) q[2];
sx q[2];
rz(-2.2714104) q[2];
sx q[2];
rz(2.0127279) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0128768) q[1];
sx q[1];
rz(-2.8059462) q[1];
sx q[1];
rz(-0.0291834) q[1];
rz(2.9167261) q[3];
sx q[3];
rz(-1.9266085) q[3];
sx q[3];
rz(-2.0840621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.31972739) q[2];
sx q[2];
rz(-1.3250019) q[2];
sx q[2];
rz(-0.74903178) q[2];
rz(-2.9876515) q[3];
sx q[3];
rz(-2.1059683) q[3];
sx q[3];
rz(-3.0416987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.446949) q[0];
sx q[0];
rz(-1.6634989) q[0];
sx q[0];
rz(-1.9248167) q[0];
rz(1.0617537) q[1];
sx q[1];
rz(-2.3371425) q[1];
sx q[1];
rz(-2.7144576) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1841284) q[0];
sx q[0];
rz(-1.9165358) q[0];
sx q[0];
rz(0.60681245) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8893625) q[2];
sx q[2];
rz(-0.67020352) q[2];
sx q[2];
rz(-0.88662749) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1041873) q[1];
sx q[1];
rz(-1.4302642) q[1];
sx q[1];
rz(-1.2804081) q[1];
rz(-pi) q[2];
rz(-1.322836) q[3];
sx q[3];
rz(-2.9313847) q[3];
sx q[3];
rz(-1.5402286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.079166807) q[2];
sx q[2];
rz(-1.7288952) q[2];
sx q[2];
rz(-1.742935) q[2];
rz(-0.22377293) q[3];
sx q[3];
rz(-0.90884915) q[3];
sx q[3];
rz(-1.5435425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1200714) q[0];
sx q[0];
rz(-1.3503617) q[0];
sx q[0];
rz(-3.0960826) q[0];
rz(-1.9728164) q[1];
sx q[1];
rz(-1.903542) q[1];
sx q[1];
rz(2.5659836) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63059124) q[0];
sx q[0];
rz(-1.5596034) q[0];
sx q[0];
rz(2.5024947) q[0];
x q[1];
rz(0.017150684) q[2];
sx q[2];
rz(-0.91020012) q[2];
sx q[2];
rz(1.2397546) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3176885) q[1];
sx q[1];
rz(-1.0270734) q[1];
sx q[1];
rz(0.25441092) q[1];
rz(-pi) q[2];
rz(3.1166385) q[3];
sx q[3];
rz(-0.91812274) q[3];
sx q[3];
rz(3.0252473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.0098972926) q[2];
sx q[2];
rz(-0.74247777) q[2];
sx q[2];
rz(-2.1843145) q[2];
rz(0.57473985) q[3];
sx q[3];
rz(-1.6114019) q[3];
sx q[3];
rz(1.3055698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15605536) q[0];
sx q[0];
rz(-0.95273459) q[0];
sx q[0];
rz(1.8804469) q[0];
rz(-0.082322923) q[1];
sx q[1];
rz(-2.0539093) q[1];
sx q[1];
rz(-1.7877158) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60786965) q[0];
sx q[0];
rz(-1.3114616) q[0];
sx q[0];
rz(1.813867) q[0];
rz(-0.48575966) q[2];
sx q[2];
rz(-1.4076715) q[2];
sx q[2];
rz(0.32147929) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8435858) q[1];
sx q[1];
rz(-1.6248967) q[1];
sx q[1];
rz(-0.49523103) q[1];
rz(-pi) q[2];
rz(2.1534377) q[3];
sx q[3];
rz(-2.372962) q[3];
sx q[3];
rz(0.60291327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7202619) q[2];
sx q[2];
rz(-2.223184) q[2];
sx q[2];
rz(1.8850231) q[2];
rz(2.4300857) q[3];
sx q[3];
rz(-1.810775) q[3];
sx q[3];
rz(0.28212696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8307777) q[0];
sx q[0];
rz(-0.84596914) q[0];
sx q[0];
rz(-2.7984483) q[0];
rz(-0.061773069) q[1];
sx q[1];
rz(-0.97007483) q[1];
sx q[1];
rz(-1.4168581) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8648711) q[0];
sx q[0];
rz(-1.7407932) q[0];
sx q[0];
rz(1.2890588) q[0];
rz(-pi) q[1];
rz(1.896865) q[2];
sx q[2];
rz(-0.89674258) q[2];
sx q[2];
rz(1.3671966) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0793905) q[1];
sx q[1];
rz(-0.50143999) q[1];
sx q[1];
rz(1.8556684) q[1];
rz(2.9067578) q[3];
sx q[3];
rz(-1.9855301) q[3];
sx q[3];
rz(-1.0339586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5938277) q[2];
sx q[2];
rz(-1.6739028) q[2];
sx q[2];
rz(1.0859547) q[2];
rz(0.91935277) q[3];
sx q[3];
rz(-1.3708401) q[3];
sx q[3];
rz(0.61387387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3980961) q[0];
sx q[0];
rz(-1.2092104) q[0];
sx q[0];
rz(-0.79291517) q[0];
rz(0.74053699) q[1];
sx q[1];
rz(-2.1426327) q[1];
sx q[1];
rz(-0.83121306) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3073458) q[0];
sx q[0];
rz(-1.1957268) q[0];
sx q[0];
rz(-2.9177279) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5969876) q[2];
sx q[2];
rz(-0.27715836) q[2];
sx q[2];
rz(-1.8089393) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.97396353) q[1];
sx q[1];
rz(-2.2669753) q[1];
sx q[1];
rz(-2.7401398) q[1];
x q[2];
rz(1.100086) q[3];
sx q[3];
rz(-1.5892972) q[3];
sx q[3];
rz(-0.61533606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.071216019) q[2];
sx q[2];
rz(-1.2717609) q[2];
sx q[2];
rz(1.1191204) q[2];
rz(1.2707155) q[3];
sx q[3];
rz(-2.0344574) q[3];
sx q[3];
rz(-1.3814111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1374461) q[0];
sx q[0];
rz(-2.9747712) q[0];
sx q[0];
rz(-1.5995837) q[0];
rz(2.1265325) q[1];
sx q[1];
rz(-1.5856182) q[1];
sx q[1];
rz(-0.17280811) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61758608) q[0];
sx q[0];
rz(-0.79389555) q[0];
sx q[0];
rz(-2.2909597) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6989231) q[2];
sx q[2];
rz(-0.73390642) q[2];
sx q[2];
rz(1.1709605) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0513251) q[1];
sx q[1];
rz(-2.6410612) q[1];
sx q[1];
rz(2.1085018) q[1];
x q[2];
rz(-0.76877131) q[3];
sx q[3];
rz(-1.1535346) q[3];
sx q[3];
rz(2.4989243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.479594) q[2];
sx q[2];
rz(-2.2981503) q[2];
sx q[2];
rz(-0.97770989) q[2];
rz(-2.5665723) q[3];
sx q[3];
rz(-1.2049371) q[3];
sx q[3];
rz(-1.2303111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32790312) q[0];
sx q[0];
rz(-0.29569018) q[0];
sx q[0];
rz(-3.1258702) q[0];
rz(-1.188259) q[1];
sx q[1];
rz(-0.63242042) q[1];
sx q[1];
rz(1.1994919) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2000113) q[0];
sx q[0];
rz(-2.0111548) q[0];
sx q[0];
rz(0.26047996) q[0];
rz(-pi) q[1];
rz(2.0616777) q[2];
sx q[2];
rz(-1.5245617) q[2];
sx q[2];
rz(0.85109988) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6024692) q[1];
sx q[1];
rz(-1.1014465) q[1];
sx q[1];
rz(-2.5825809) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1114798) q[3];
sx q[3];
rz(-0.62255854) q[3];
sx q[3];
rz(-0.42296577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1477995) q[2];
sx q[2];
rz(-1.9027998) q[2];
sx q[2];
rz(-1.3851059) q[2];
rz(-0.6238474) q[3];
sx q[3];
rz(-1.8892989) q[3];
sx q[3];
rz(1.0998211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95595908) q[0];
sx q[0];
rz(-2.3210242) q[0];
sx q[0];
rz(-2.4483335) q[0];
rz(0.26501003) q[1];
sx q[1];
rz(-0.82844228) q[1];
sx q[1];
rz(-1.5230491) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.023848195) q[0];
sx q[0];
rz(-1.5597514) q[0];
sx q[0];
rz(-2.3680229) q[0];
rz(-pi) q[1];
rz(-0.9829282) q[2];
sx q[2];
rz(-2.1097906) q[2];
sx q[2];
rz(1.3608152) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0761145) q[1];
sx q[1];
rz(-1.5947184) q[1];
sx q[1];
rz(0.33965276) q[1];
rz(-pi) q[2];
rz(1.4632439) q[3];
sx q[3];
rz(-0.58379025) q[3];
sx q[3];
rz(1.8369758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.30274621) q[2];
sx q[2];
rz(-1.5134209) q[2];
sx q[2];
rz(-1.6839074) q[2];
rz(0.95528209) q[3];
sx q[3];
rz(-0.49946076) q[3];
sx q[3];
rz(0.83520755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7569698) q[0];
sx q[0];
rz(-0.56461016) q[0];
sx q[0];
rz(-0.48267522) q[0];
rz(0.92974281) q[1];
sx q[1];
rz(-1.0044121) q[1];
sx q[1];
rz(0.39628705) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1866465) q[0];
sx q[0];
rz(-2.5781693) q[0];
sx q[0];
rz(2.4231829) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9626289) q[2];
sx q[2];
rz(-1.4216627) q[2];
sx q[2];
rz(2.3863132) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0490108) q[1];
sx q[1];
rz(-1.7515469) q[1];
sx q[1];
rz(0.84047079) q[1];
rz(-pi) q[2];
rz(0.72209218) q[3];
sx q[3];
rz(-0.51755899) q[3];
sx q[3];
rz(-0.90921569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.79260176) q[2];
sx q[2];
rz(-1.7020117) q[2];
sx q[2];
rz(-0.3271884) q[2];
rz(2.3606825) q[3];
sx q[3];
rz(-1.9802997) q[3];
sx q[3];
rz(-1.4769295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1031716) q[0];
sx q[0];
rz(-2.1506943) q[0];
sx q[0];
rz(-2.8839169) q[0];
rz(-2.7720263) q[1];
sx q[1];
rz(-2.259544) q[1];
sx q[1];
rz(2.4824711) q[1];
rz(-2.4482881) q[2];
sx q[2];
rz(-1.605576) q[2];
sx q[2];
rz(-1.3645542) q[2];
rz(2.1791069) q[3];
sx q[3];
rz(-1.8909834) q[3];
sx q[3];
rz(-3.0743619) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
