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
rz(-1.3136366) q[0];
sx q[0];
rz(-3.0744636) q[0];
sx q[0];
rz(-0.022775291) q[0];
rz(1.0442806) q[1];
sx q[1];
rz(-2.2106946) q[1];
sx q[1];
rz(-3.1276303) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9305796) q[0];
sx q[0];
rz(-1.6114116) q[0];
sx q[0];
rz(-2.3848349) q[0];
rz(-pi) q[1];
rz(0.44273744) q[2];
sx q[2];
rz(-1.95089) q[2];
sx q[2];
rz(0.9600823) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3252445) q[1];
sx q[1];
rz(-0.26953408) q[1];
sx q[1];
rz(2.4052038) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9223677) q[3];
sx q[3];
rz(-1.1446125) q[3];
sx q[3];
rz(3.1015729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.969101) q[2];
sx q[2];
rz(-1.3250985) q[2];
sx q[2];
rz(1.8005523) q[2];
rz(-1.7155044) q[3];
sx q[3];
rz(-1.1172349) q[3];
sx q[3];
rz(-2.8360227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39389998) q[0];
sx q[0];
rz(-1.9992398) q[0];
sx q[0];
rz(-0.72545141) q[0];
rz(-0.91616383) q[1];
sx q[1];
rz(-2.3326645) q[1];
sx q[1];
rz(-2.5221672) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40901285) q[0];
sx q[0];
rz(-2.469099) q[0];
sx q[0];
rz(-0.88100638) q[0];
rz(-2.5079501) q[2];
sx q[2];
rz(-0.85553193) q[2];
sx q[2];
rz(1.9618386) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.933085) q[1];
sx q[1];
rz(-1.2080482) q[1];
sx q[1];
rz(-2.0278992) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1257915) q[3];
sx q[3];
rz(-1.8280622) q[3];
sx q[3];
rz(2.1580838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.82975036) q[2];
sx q[2];
rz(-2.1205015) q[2];
sx q[2];
rz(-0.96168438) q[2];
rz(-2.035615) q[3];
sx q[3];
rz(-1.032369) q[3];
sx q[3];
rz(0.34524125) q[3];
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
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2894534) q[0];
sx q[0];
rz(-1.1281321) q[0];
sx q[0];
rz(-1.7732675) q[0];
rz(0.013956919) q[1];
sx q[1];
rz(-2.0266666) q[1];
sx q[1];
rz(-1.7272635) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089979261) q[0];
sx q[0];
rz(-0.3246626) q[0];
sx q[0];
rz(-1.8568296) q[0];
rz(0.0711381) q[2];
sx q[2];
rz(-1.2025598) q[2];
sx q[2];
rz(-0.28063831) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5217375) q[1];
sx q[1];
rz(-0.69112366) q[1];
sx q[1];
rz(-2.28165) q[1];
rz(-pi) q[2];
rz(-0.39140626) q[3];
sx q[3];
rz(-1.3308755) q[3];
sx q[3];
rz(2.4521884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3783375) q[2];
sx q[2];
rz(-1.1474643) q[2];
sx q[2];
rz(-2.9249127) q[2];
rz(-2.2583151) q[3];
sx q[3];
rz(-2.7933385) q[3];
sx q[3];
rz(-1.1368375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9388409) q[0];
sx q[0];
rz(-1.0067679) q[0];
sx q[0];
rz(2.3034565) q[0];
rz(-0.54191598) q[1];
sx q[1];
rz(-0.37626615) q[1];
sx q[1];
rz(0.045104973) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7823001) q[0];
sx q[0];
rz(-1.4385537) q[0];
sx q[0];
rz(0.6036997) q[0];
x q[1];
rz(-2.0531282) q[2];
sx q[2];
rz(-1.3399897) q[2];
sx q[2];
rz(-2.7969691) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.038626898) q[1];
sx q[1];
rz(-0.47768394) q[1];
sx q[1];
rz(-1.4250907) q[1];
x q[2];
rz(1.6691505) q[3];
sx q[3];
rz(-1.0788267) q[3];
sx q[3];
rz(1.4592001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1230459) q[2];
sx q[2];
rz(-0.16877731) q[2];
sx q[2];
rz(-2.4875842) q[2];
rz(2.5096014) q[3];
sx q[3];
rz(-1.4258823) q[3];
sx q[3];
rz(-2.8051207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80394799) q[0];
sx q[0];
rz(-2.4027282) q[0];
sx q[0];
rz(-1.3496572) q[0];
rz(-2.3140287) q[1];
sx q[1];
rz(-2.5065828) q[1];
sx q[1];
rz(2.6571224) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.890936) q[0];
sx q[0];
rz(-1.2095426) q[0];
sx q[0];
rz(2.7879232) q[0];
rz(-pi) q[1];
rz(-3.0516009) q[2];
sx q[2];
rz(-1.0183739) q[2];
sx q[2];
rz(-1.8664139) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3282991) q[1];
sx q[1];
rz(-0.44078953) q[1];
sx q[1];
rz(-0.30118024) q[1];
rz(-0.46273939) q[3];
sx q[3];
rz(-2.2125508) q[3];
sx q[3];
rz(-0.79824191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1822002) q[2];
sx q[2];
rz(-1.0687989) q[2];
sx q[2];
rz(2.5051795) q[2];
rz(-0.098879769) q[3];
sx q[3];
rz(-1.585588) q[3];
sx q[3];
rz(1.4720565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32475489) q[0];
sx q[0];
rz(-0.88961283) q[0];
sx q[0];
rz(-1.1725934) q[0];
rz(1.532754) q[1];
sx q[1];
rz(-1.0120665) q[1];
sx q[1];
rz(-0.00072678725) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1689735) q[0];
sx q[0];
rz(-1.0797636) q[0];
sx q[0];
rz(-2.6537077) q[0];
rz(1.4719719) q[2];
sx q[2];
rz(-1.3512058) q[2];
sx q[2];
rz(-1.7864986) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.24405133) q[1];
sx q[1];
rz(-1.6836299) q[1];
sx q[1];
rz(2.807602) q[1];
rz(-0.57925333) q[3];
sx q[3];
rz(-1.4046852) q[3];
sx q[3];
rz(0.36234713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8751004) q[2];
sx q[2];
rz(-0.88783395) q[2];
sx q[2];
rz(2.9317741) q[2];
rz(0.80875129) q[3];
sx q[3];
rz(-2.5918312) q[3];
sx q[3];
rz(-0.89890283) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4631735) q[0];
sx q[0];
rz(-1.9371978) q[0];
sx q[0];
rz(3.0294982) q[0];
rz(0.049292715) q[1];
sx q[1];
rz(-0.53105989) q[1];
sx q[1];
rz(1.3370399) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61510117) q[0];
sx q[0];
rz(-1.3630629) q[0];
sx q[0];
rz(-0.92089842) q[0];
x q[1];
rz(-2.078901) q[2];
sx q[2];
rz(-1.6828949) q[2];
sx q[2];
rz(1.0401806) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.783205) q[1];
sx q[1];
rz(-1.9175944) q[1];
sx q[1];
rz(1.5601182) q[1];
rz(-pi) q[2];
rz(1.8932166) q[3];
sx q[3];
rz(-2.5395576) q[3];
sx q[3];
rz(0.69794929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2992531) q[2];
sx q[2];
rz(-2.3167593) q[2];
sx q[2];
rz(0.71005589) q[2];
rz(0.0013141343) q[3];
sx q[3];
rz(-2.240286) q[3];
sx q[3];
rz(-0.60091758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8173219) q[0];
sx q[0];
rz(-1.7621499) q[0];
sx q[0];
rz(-2.5589909) q[0];
rz(-0.6800037) q[1];
sx q[1];
rz(-0.99907196) q[1];
sx q[1];
rz(1.2635788) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9506142) q[0];
sx q[0];
rz(-0.54784566) q[0];
sx q[0];
rz(0.56444278) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.81496691) q[2];
sx q[2];
rz(-0.294058) q[2];
sx q[2];
rz(0.61758274) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3547825) q[1];
sx q[1];
rz(-2.218333) q[1];
sx q[1];
rz(-1.3446273) q[1];
rz(-0.42383343) q[3];
sx q[3];
rz(-1.5409711) q[3];
sx q[3];
rz(0.78889293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6363643) q[2];
sx q[2];
rz(-1.4437081) q[2];
sx q[2];
rz(-0.12602028) q[2];
rz(-1.5382918) q[3];
sx q[3];
rz(-2.3647629) q[3];
sx q[3];
rz(-2.7974424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4890323) q[0];
sx q[0];
rz(-1.6283988) q[0];
sx q[0];
rz(-1.162758) q[0];
rz(-1.8310422) q[1];
sx q[1];
rz(-2.4030011) q[1];
sx q[1];
rz(1.383925) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21400586) q[0];
sx q[0];
rz(-1.5259552) q[0];
sx q[0];
rz(-1.2311185) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6430364) q[2];
sx q[2];
rz(-1.6933228) q[2];
sx q[2];
rz(0.83736698) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0392929) q[1];
sx q[1];
rz(-0.96233923) q[1];
sx q[1];
rz(-0.14157544) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8683327) q[3];
sx q[3];
rz(-0.14500824) q[3];
sx q[3];
rz(2.6535237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3414574) q[2];
sx q[2];
rz(-1.2722445) q[2];
sx q[2];
rz(0.47194353) q[2];
rz(2.3927472) q[3];
sx q[3];
rz(-2.2180836) q[3];
sx q[3];
rz(2.3239465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6227459) q[0];
sx q[0];
rz(-2.2846344) q[0];
sx q[0];
rz(2.6527606) q[0];
rz(-0.34293109) q[1];
sx q[1];
rz(-0.88645005) q[1];
sx q[1];
rz(-2.0874646) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80741548) q[0];
sx q[0];
rz(-1.6082967) q[0];
sx q[0];
rz(1.6223909) q[0];
rz(2.6661886) q[2];
sx q[2];
rz(-0.66988551) q[2];
sx q[2];
rz(0.49150447) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3588402) q[1];
sx q[1];
rz(-0.61780518) q[1];
sx q[1];
rz(-2.0899523) q[1];
x q[2];
rz(-1.6659237) q[3];
sx q[3];
rz(-0.47387487) q[3];
sx q[3];
rz(-2.9364862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1663345) q[2];
sx q[2];
rz(-2.9322093) q[2];
sx q[2];
rz(-0.091910467) q[2];
rz(2.6722243) q[3];
sx q[3];
rz(-2.2769603) q[3];
sx q[3];
rz(-0.11274591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54409201) q[0];
sx q[0];
rz(-1.8202029) q[0];
sx q[0];
rz(1.2910917) q[0];
rz(2.3757833) q[1];
sx q[1];
rz(-1.7821093) q[1];
sx q[1];
rz(-1.2761188) q[1];
rz(0.2551887) q[2];
sx q[2];
rz(-1.219426) q[2];
sx q[2];
rz(0.30795369) q[2];
rz(0.53743955) q[3];
sx q[3];
rz(-1.3466571) q[3];
sx q[3];
rz(1.3206645) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
