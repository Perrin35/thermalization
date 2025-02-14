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
rz(2.9003484) q[0];
sx q[0];
rz(3.3527346) q[0];
sx q[0];
rz(9.0169173) q[0];
rz(-1.7653699) q[1];
sx q[1];
rz(5.0566109) q[1];
sx q[1];
rz(15.772718) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0000251) q[0];
sx q[0];
rz(-1.5760826) q[0];
sx q[0];
rz(-1.686211) q[0];
rz(-0.32081713) q[2];
sx q[2];
rz(-2.1324369) q[2];
sx q[2];
rz(-0.80036847) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3169466) q[1];
sx q[1];
rz(-1.4295409) q[1];
sx q[1];
rz(1.5596175) q[1];
rz(-pi) q[2];
x q[2];
rz(1.466732) q[3];
sx q[3];
rz(-0.95115653) q[3];
sx q[3];
rz(2.5980169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.18225081) q[2];
sx q[2];
rz(-2.0472517) q[2];
sx q[2];
rz(-1.135929) q[2];
rz(0.7302537) q[3];
sx q[3];
rz(-1.3948995) q[3];
sx q[3];
rz(1.5203389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.092904329) q[0];
sx q[0];
rz(-0.95656675) q[0];
sx q[0];
rz(3.0862869) q[0];
rz(2.7703908) q[1];
sx q[1];
rz(-0.35938811) q[1];
sx q[1];
rz(0.41195437) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7474351) q[0];
sx q[0];
rz(-2.0965946) q[0];
sx q[0];
rz(1.6461685) q[0];
rz(-pi) q[1];
rz(2.2949176) q[2];
sx q[2];
rz(-2.7077423) q[2];
sx q[2];
rz(2.2920319) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4476599) q[1];
sx q[1];
rz(-2.3946163) q[1];
sx q[1];
rz(-0.20836094) q[1];
x q[2];
rz(1.0874416) q[3];
sx q[3];
rz(-2.4341741) q[3];
sx q[3];
rz(-0.055295769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.61341316) q[2];
sx q[2];
rz(-1.6772905) q[2];
sx q[2];
rz(0.8826274) q[2];
rz(1.4334076) q[3];
sx q[3];
rz(-1.2127168) q[3];
sx q[3];
rz(2.9679756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16647896) q[0];
sx q[0];
rz(-3.1298895) q[0];
sx q[0];
rz(0.56667462) q[0];
rz(0.4176248) q[1];
sx q[1];
rz(-1.6028701) q[1];
sx q[1];
rz(-1.3272939) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9629563) q[0];
sx q[0];
rz(-1.6554852) q[0];
sx q[0];
rz(-0.87015193) q[0];
rz(-1.3603052) q[2];
sx q[2];
rz(-2.5062525) q[2];
sx q[2];
rz(0.47266211) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7975038) q[1];
sx q[1];
rz(-2.833539) q[1];
sx q[1];
rz(-1.5750242) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0762844) q[3];
sx q[3];
rz(-2.3479241) q[3];
sx q[3];
rz(-0.46773673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7495482) q[2];
sx q[2];
rz(-1.0907402) q[2];
sx q[2];
rz(-0.90847477) q[2];
rz(-2.3467017) q[3];
sx q[3];
rz(-1.1029693) q[3];
sx q[3];
rz(0.046549646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6598776) q[0];
sx q[0];
rz(-2.2292723) q[0];
sx q[0];
rz(0.6657486) q[0];
rz(2.2944229) q[1];
sx q[1];
rz(-0.24996346) q[1];
sx q[1];
rz(2.7319103) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6981229) q[0];
sx q[0];
rz(-1.4181678) q[0];
sx q[0];
rz(0.49907719) q[0];
rz(-pi) q[1];
rz(-0.3607765) q[2];
sx q[2];
rz(-2.2998126) q[2];
sx q[2];
rz(-1.3363802) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6050501) q[1];
sx q[1];
rz(-0.66075051) q[1];
sx q[1];
rz(-1.9416757) q[1];
rz(3.112875) q[3];
sx q[3];
rz(-2.014694) q[3];
sx q[3];
rz(-0.2335399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.72982558) q[2];
sx q[2];
rz(-2.4335786) q[2];
sx q[2];
rz(1.3789619) q[2];
rz(1.1369368) q[3];
sx q[3];
rz(-1.3552908) q[3];
sx q[3];
rz(0.8374477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9627422) q[0];
sx q[0];
rz(-2.0806291) q[0];
sx q[0];
rz(0.41819292) q[0];
rz(-1.7182982) q[1];
sx q[1];
rz(-1.1885252) q[1];
sx q[1];
rz(-1.4994015) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2432775) q[0];
sx q[0];
rz(-2.0137798) q[0];
sx q[0];
rz(-1.8421296) q[0];
rz(-pi) q[1];
rz(-2.8078305) q[2];
sx q[2];
rz(-1.2818206) q[2];
sx q[2];
rz(-0.0043650345) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.114257) q[1];
sx q[1];
rz(-1.789195) q[1];
sx q[1];
rz(-0.059192358) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0357596) q[3];
sx q[3];
rz(-0.14497193) q[3];
sx q[3];
rz(1.5752058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8666009) q[2];
sx q[2];
rz(-1.9860705) q[2];
sx q[2];
rz(0.17821136) q[2];
rz(-0.37402672) q[3];
sx q[3];
rz(-2.1686797) q[3];
sx q[3];
rz(-1.9988352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6571534) q[0];
sx q[0];
rz(-1.8171808) q[0];
sx q[0];
rz(-1.665218) q[0];
rz(0.58552512) q[1];
sx q[1];
rz(-0.89591566) q[1];
sx q[1];
rz(0.73385986) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4910262) q[0];
sx q[0];
rz(-2.0163476) q[0];
sx q[0];
rz(-0.70029052) q[0];
x q[1];
rz(-2.8346905) q[2];
sx q[2];
rz(-1.2595121) q[2];
sx q[2];
rz(0.5615304) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.1230525) q[1];
sx q[1];
rz(-0.55852671) q[1];
sx q[1];
rz(-1.194242) q[1];
rz(-pi) q[2];
rz(2.0229983) q[3];
sx q[3];
rz(-1.3824995) q[3];
sx q[3];
rz(1.2533497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0840941) q[2];
sx q[2];
rz(-1.7666662) q[2];
sx q[2];
rz(1.499184) q[2];
rz(-1.621014) q[3];
sx q[3];
rz(-2.4937544) q[3];
sx q[3];
rz(0.15629855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75477377) q[0];
sx q[0];
rz(-0.91803011) q[0];
sx q[0];
rz(1.7479489) q[0];
rz(-0.58760324) q[1];
sx q[1];
rz(-1.6410476) q[1];
sx q[1];
rz(-2.844152) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35703941) q[0];
sx q[0];
rz(-1.2149356) q[0];
sx q[0];
rz(0.83609348) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1226899) q[2];
sx q[2];
rz(-1.3002965) q[2];
sx q[2];
rz(0.30735415) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6828281) q[1];
sx q[1];
rz(-1.5804844) q[1];
sx q[1];
rz(-2.2504456) q[1];
x q[2];
rz(0.92771156) q[3];
sx q[3];
rz(-1.3731908) q[3];
sx q[3];
rz(-1.2090982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8224767) q[2];
sx q[2];
rz(-1.3603223) q[2];
sx q[2];
rz(-0.3817257) q[2];
rz(0.60025674) q[3];
sx q[3];
rz(-2.468942) q[3];
sx q[3];
rz(-2.8762347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47939077) q[0];
sx q[0];
rz(-1.6265765) q[0];
sx q[0];
rz(1.3375244) q[0];
rz(3.1321101) q[1];
sx q[1];
rz(-0.9597221) q[1];
sx q[1];
rz(-1.3710075) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.079938) q[0];
sx q[0];
rz(-1.2324872) q[0];
sx q[0];
rz(-1.4251955) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6422136) q[2];
sx q[2];
rz(-2.1771095) q[2];
sx q[2];
rz(-2.2057836) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7035455) q[1];
sx q[1];
rz(-0.79709541) q[1];
sx q[1];
rz(0.17532562) q[1];
rz(1.6608606) q[3];
sx q[3];
rz(-1.1582631) q[3];
sx q[3];
rz(-1.5145258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3235772) q[2];
sx q[2];
rz(-0.90138268) q[2];
sx q[2];
rz(-2.9991007) q[2];
rz(-0.88322181) q[3];
sx q[3];
rz(-1.3775237) q[3];
sx q[3];
rz(-1.4628598) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0505117) q[0];
sx q[0];
rz(-0.091878042) q[0];
sx q[0];
rz(-1.8930513) q[0];
rz(2.7342791) q[1];
sx q[1];
rz(-1.3469603) q[1];
sx q[1];
rz(-1.9479082) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4536205) q[0];
sx q[0];
rz(-1.2415452) q[0];
sx q[0];
rz(-0.27197522) q[0];
x q[1];
rz(-2.418708) q[2];
sx q[2];
rz(-2.2396846) q[2];
sx q[2];
rz(-1.0793874) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9858169) q[1];
sx q[1];
rz(-1.078842) q[1];
sx q[1];
rz(1.953275) q[1];
x q[2];
rz(-1.3109929) q[3];
sx q[3];
rz(-1.570829) q[3];
sx q[3];
rz(0.81476975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.074389) q[2];
sx q[2];
rz(-0.95778242) q[2];
sx q[2];
rz(-2.0286782) q[2];
rz(-3.0985966) q[3];
sx q[3];
rz(-1.6578511) q[3];
sx q[3];
rz(-2.903741) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.064903108) q[0];
sx q[0];
rz(-1.4142798) q[0];
sx q[0];
rz(-2.9458556) q[0];
rz(-1.2204569) q[1];
sx q[1];
rz(-0.38206044) q[1];
sx q[1];
rz(2.1641796) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.138151) q[0];
sx q[0];
rz(-1.312698) q[0];
sx q[0];
rz(-1.6736223) q[0];
x q[1];
rz(1.8711817) q[2];
sx q[2];
rz(-1.7160048) q[2];
sx q[2];
rz(1.6329488) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7070476) q[1];
sx q[1];
rz(-2.0710398) q[1];
sx q[1];
rz(-2.2342199) q[1];
rz(-pi) q[2];
rz(-1.8582088) q[3];
sx q[3];
rz(-1.948878) q[3];
sx q[3];
rz(-1.9668129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0105373) q[2];
sx q[2];
rz(-0.26630339) q[2];
sx q[2];
rz(1.4492501) q[2];
rz(0.87861711) q[3];
sx q[3];
rz(-1.0613469) q[3];
sx q[3];
rz(-1.3531468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5129678) q[0];
sx q[0];
rz(-1.9504564) q[0];
sx q[0];
rz(1.7497028) q[0];
rz(1.7616918) q[1];
sx q[1];
rz(-1.5329783) q[1];
sx q[1];
rz(3.0402532) q[1];
rz(-1.7529132) q[2];
sx q[2];
rz(-1.7487329) q[2];
sx q[2];
rz(2.6404811) q[2];
rz(2.8158549) q[3];
sx q[3];
rz(-0.67254638) q[3];
sx q[3];
rz(-1.2253958) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
