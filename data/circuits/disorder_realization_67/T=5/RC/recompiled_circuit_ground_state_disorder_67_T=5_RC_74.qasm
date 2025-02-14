OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9031653) q[0];
sx q[0];
rz(2.981346) q[0];
sx q[0];
rz(7.5709406) q[0];
rz(2.2924478) q[1];
sx q[1];
rz(-1.949911) q[1];
sx q[1];
rz(2.3828659) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.01837) q[0];
sx q[0];
rz(-1.8022853) q[0];
sx q[0];
rz(2.8442848) q[0];
x q[1];
rz(2.0082194) q[2];
sx q[2];
rz(-2.6133399) q[2];
sx q[2];
rz(-0.6483486) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.78658653) q[1];
sx q[1];
rz(-2.2753582) q[1];
sx q[1];
rz(-0.48272065) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.82697273) q[3];
sx q[3];
rz(-1.0543514) q[3];
sx q[3];
rz(1.6940821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8410926) q[2];
sx q[2];
rz(-0.79885834) q[2];
sx q[2];
rz(-2.63499) q[2];
rz(1.6182342) q[3];
sx q[3];
rz(-0.62464276) q[3];
sx q[3];
rz(3.0860331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.8697206) q[0];
sx q[0];
rz(-2.6756918) q[0];
sx q[0];
rz(-3.077935) q[0];
rz(1.9438538) q[1];
sx q[1];
rz(-1.627219) q[1];
sx q[1];
rz(2.1228085) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1757024) q[0];
sx q[0];
rz(-1.4251484) q[0];
sx q[0];
rz(-0.19840981) q[0];
rz(-1.8070142) q[2];
sx q[2];
rz(-1.902632) q[2];
sx q[2];
rz(-2.163909) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0662811) q[1];
sx q[1];
rz(-2.3763658) q[1];
sx q[1];
rz(-0.40513961) q[1];
rz(-pi) q[2];
rz(-0.81786675) q[3];
sx q[3];
rz(-1.6685155) q[3];
sx q[3];
rz(2.1678501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7055052) q[2];
sx q[2];
rz(-2.5477396) q[2];
sx q[2];
rz(0.84211055) q[2];
rz(2.0841133) q[3];
sx q[3];
rz(-0.92856854) q[3];
sx q[3];
rz(-1.1697945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24404003) q[0];
sx q[0];
rz(-1.0201447) q[0];
sx q[0];
rz(-1.1919588) q[0];
rz(2.2167218) q[1];
sx q[1];
rz(-0.7083188) q[1];
sx q[1];
rz(-1.8398197) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.155543) q[0];
sx q[0];
rz(-1.0240004) q[0];
sx q[0];
rz(0.56122551) q[0];
rz(-pi) q[1];
rz(-1.6100175) q[2];
sx q[2];
rz(-1.3717781) q[2];
sx q[2];
rz(1.4250371) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8346602) q[1];
sx q[1];
rz(-1.3921157) q[1];
sx q[1];
rz(2.6684922) q[1];
rz(-pi) q[2];
rz(0.93940063) q[3];
sx q[3];
rz(-2.7223848) q[3];
sx q[3];
rz(-1.4155751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3785582) q[2];
sx q[2];
rz(-3.0578461) q[2];
sx q[2];
rz(-0.087510022) q[2];
rz(-1.3148426) q[3];
sx q[3];
rz(-1.0564691) q[3];
sx q[3];
rz(2.1480613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73115504) q[0];
sx q[0];
rz(-1.9696099) q[0];
sx q[0];
rz(0.80899578) q[0];
rz(2.1058829) q[1];
sx q[1];
rz(-2.6782942) q[1];
sx q[1];
rz(-1.0332003) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22543487) q[0];
sx q[0];
rz(-1.5317129) q[0];
sx q[0];
rz(2.2351511) q[0];
rz(-2.2351) q[2];
sx q[2];
rz(-2.4943753) q[2];
sx q[2];
rz(-2.3075019) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.82419187) q[1];
sx q[1];
rz(-2.6397974) q[1];
sx q[1];
rz(0.3806033) q[1];
rz(2.4304588) q[3];
sx q[3];
rz(-0.67340771) q[3];
sx q[3];
rz(-2.8720958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0741299) q[2];
sx q[2];
rz(-0.16572696) q[2];
sx q[2];
rz(1.5376252) q[2];
rz(1.4175203) q[3];
sx q[3];
rz(-1.1105024) q[3];
sx q[3];
rz(-2.0541151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9017225) q[0];
sx q[0];
rz(-0.19342315) q[0];
sx q[0];
rz(-1.9523917) q[0];
rz(-0.72422782) q[1];
sx q[1];
rz(-1.8502356) q[1];
sx q[1];
rz(1.474818) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14872197) q[0];
sx q[0];
rz(-1.9970511) q[0];
sx q[0];
rz(0.54963407) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3442845) q[2];
sx q[2];
rz(-1.6835137) q[2];
sx q[2];
rz(-1.1657451) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.96383038) q[1];
sx q[1];
rz(-0.46198598) q[1];
sx q[1];
rz(2.3038008) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.83051564) q[3];
sx q[3];
rz(-0.6336824) q[3];
sx q[3];
rz(-0.049402852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7597947) q[2];
sx q[2];
rz(-0.90961027) q[2];
sx q[2];
rz(2.6749715) q[2];
rz(-1.4383379) q[3];
sx q[3];
rz(-1.8432901) q[3];
sx q[3];
rz(-0.94021016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71317116) q[0];
sx q[0];
rz(-2.3155825) q[0];
sx q[0];
rz(3.1312422) q[0];
rz(1.5230644) q[1];
sx q[1];
rz(-1.7489988) q[1];
sx q[1];
rz(1.7656322) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3165163) q[0];
sx q[0];
rz(-2.7211271) q[0];
sx q[0];
rz(-2.1964873) q[0];
rz(1.7014241) q[2];
sx q[2];
rz(-1.3816026) q[2];
sx q[2];
rz(1.4551324) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.048307) q[1];
sx q[1];
rz(-1.7518338) q[1];
sx q[1];
rz(-0.50706086) q[1];
rz(-pi) q[2];
x q[2];
rz(0.74919219) q[3];
sx q[3];
rz(-1.0674879) q[3];
sx q[3];
rz(2.1319413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4202262) q[2];
sx q[2];
rz(-1.0488291) q[2];
sx q[2];
rz(-0.72727195) q[2];
rz(2.6266802) q[3];
sx q[3];
rz(-1.3313096) q[3];
sx q[3];
rz(1.4179199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10776831) q[0];
sx q[0];
rz(-1.6078147) q[0];
sx q[0];
rz(-2.7346101) q[0];
rz(-0.96626967) q[1];
sx q[1];
rz(-2.2844908) q[1];
sx q[1];
rz(-0.13872096) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16628597) q[0];
sx q[0];
rz(-1.5353256) q[0];
sx q[0];
rz(-1.6871638) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3499603) q[2];
sx q[2];
rz(-1.8419617) q[2];
sx q[2];
rz(1.3214932) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3073716) q[1];
sx q[1];
rz(-1.929793) q[1];
sx q[1];
rz(-1.2008083) q[1];
rz(-2.9937708) q[3];
sx q[3];
rz(-1.9114219) q[3];
sx q[3];
rz(0.56440777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8653284) q[2];
sx q[2];
rz(-1.7636969) q[2];
sx q[2];
rz(2.830937) q[2];
rz(-1.3079414) q[3];
sx q[3];
rz(-1.5111978) q[3];
sx q[3];
rz(-2.7697897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0528316) q[0];
sx q[0];
rz(-2.796687) q[0];
sx q[0];
rz(-3.0533691) q[0];
rz(3.131033) q[1];
sx q[1];
rz(-1.5225531) q[1];
sx q[1];
rz(0.47058502) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39401606) q[0];
sx q[0];
rz(-2.2758099) q[0];
sx q[0];
rz(0.88253077) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9520985) q[2];
sx q[2];
rz(-2.3255286) q[2];
sx q[2];
rz(0.59743728) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7631665) q[1];
sx q[1];
rz(-1.4861123) q[1];
sx q[1];
rz(1.8316818) q[1];
x q[2];
rz(1.2180738) q[3];
sx q[3];
rz(-1.439487) q[3];
sx q[3];
rz(-2.2985947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.18870246) q[2];
sx q[2];
rz(-0.95658335) q[2];
sx q[2];
rz(-1.0799705) q[2];
rz(-2.2522669) q[3];
sx q[3];
rz(-1.5998799) q[3];
sx q[3];
rz(1.5842452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(2.476986) q[0];
sx q[0];
rz(-0.38170686) q[0];
sx q[0];
rz(0.81740022) q[0];
rz(0.74642247) q[1];
sx q[1];
rz(-0.67459977) q[1];
sx q[1];
rz(0.93961632) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22720756) q[0];
sx q[0];
rz(-0.017134754) q[0];
sx q[0];
rz(-1.5428154) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4822761) q[2];
sx q[2];
rz(-0.43956471) q[2];
sx q[2];
rz(2.2826113) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6691362) q[1];
sx q[1];
rz(-1.3208773) q[1];
sx q[1];
rz(-0.12052287) q[1];
x q[2];
rz(2.0321531) q[3];
sx q[3];
rz(-0.45518866) q[3];
sx q[3];
rz(-0.92092848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8115936) q[2];
sx q[2];
rz(-2.1758695) q[2];
sx q[2];
rz(-2.7962371) q[2];
rz(1.5251478) q[3];
sx q[3];
rz(-1.350178) q[3];
sx q[3];
rz(-2.7143872) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78373194) q[0];
sx q[0];
rz(-1.1564199) q[0];
sx q[0];
rz(0.082948908) q[0];
rz(0.16231617) q[1];
sx q[1];
rz(-1.4444618) q[1];
sx q[1];
rz(0.69581318) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0950553) q[0];
sx q[0];
rz(-0.2827321) q[0];
sx q[0];
rz(-0.15034349) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1386069) q[2];
sx q[2];
rz(-1.4415359) q[2];
sx q[2];
rz(-2.7832727) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.79244991) q[1];
sx q[1];
rz(-2.6969243) q[1];
sx q[1];
rz(-0.044388958) q[1];
rz(-pi) q[2];
x q[2];
rz(0.331075) q[3];
sx q[3];
rz(-0.73270117) q[3];
sx q[3];
rz(-0.98950451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.42056981) q[2];
sx q[2];
rz(-0.15779725) q[2];
sx q[2];
rz(1.921462) q[2];
rz(-0.57341352) q[3];
sx q[3];
rz(-1.8107332) q[3];
sx q[3];
rz(-2.8964608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86380105) q[0];
sx q[0];
rz(-1.6391123) q[0];
sx q[0];
rz(3.0085051) q[0];
rz(-3.1387023) q[1];
sx q[1];
rz(-0.4160226) q[1];
sx q[1];
rz(-0.70152534) q[1];
rz(1.3221424) q[2];
sx q[2];
rz(-1.7380309) q[2];
sx q[2];
rz(-2.0792014) q[2];
rz(3.0843432) q[3];
sx q[3];
rz(-2.3404239) q[3];
sx q[3];
rz(-0.41062582) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
