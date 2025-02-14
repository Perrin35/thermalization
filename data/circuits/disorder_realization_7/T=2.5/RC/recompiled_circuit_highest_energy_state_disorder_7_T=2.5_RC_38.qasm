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
rz(-0.67686358) q[0];
sx q[0];
rz(-1.6114177) q[0];
sx q[0];
rz(-2.8443008) q[0];
rz(0.77904207) q[1];
sx q[1];
rz(-0.160633) q[1];
sx q[1];
rz(1.4359441) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6061062) q[0];
sx q[0];
rz(-2.1155042) q[0];
sx q[0];
rz(0.60118712) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0129635) q[2];
sx q[2];
rz(-1.8278215) q[2];
sx q[2];
rz(2.2976053) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2816726) q[1];
sx q[1];
rz(-1.249673) q[1];
sx q[1];
rz(-0.74076498) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7706031) q[3];
sx q[3];
rz(-1.8733896) q[3];
sx q[3];
rz(1.2452919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5433189) q[2];
sx q[2];
rz(-0.42793772) q[2];
sx q[2];
rz(-0.16779009) q[2];
rz(-2.0134036) q[3];
sx q[3];
rz(-1.5690683) q[3];
sx q[3];
rz(-1.9523841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.096864916) q[0];
sx q[0];
rz(-2.4714578) q[0];
sx q[0];
rz(1.0750394) q[0];
rz(1.7754405) q[1];
sx q[1];
rz(-1.3864044) q[1];
sx q[1];
rz(1.3137821) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15720651) q[0];
sx q[0];
rz(-1.3493378) q[0];
sx q[0];
rz(-0.33672543) q[0];
rz(-pi) q[1];
rz(-0.18630917) q[2];
sx q[2];
rz(-1.7373475) q[2];
sx q[2];
rz(-0.09504091) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.35006501) q[1];
sx q[1];
rz(-0.72775562) q[1];
sx q[1];
rz(1.6863053) q[1];
rz(1.2398534) q[3];
sx q[3];
rz(-1.6972491) q[3];
sx q[3];
rz(3.0883873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4904867) q[2];
sx q[2];
rz(-1.7650812) q[2];
sx q[2];
rz(-0.42438486) q[2];
rz(0.73244798) q[3];
sx q[3];
rz(-1.6583574) q[3];
sx q[3];
rz(-1.463416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3378147) q[0];
sx q[0];
rz(-2.4265899) q[0];
sx q[0];
rz(-2.6388229) q[0];
rz(-2.1979507) q[1];
sx q[1];
rz(-0.6424526) q[1];
sx q[1];
rz(2.3263993) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73683263) q[0];
sx q[0];
rz(-1.6175999) q[0];
sx q[0];
rz(-0.18808774) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.11932245) q[2];
sx q[2];
rz(-1.4652958) q[2];
sx q[2];
rz(-1.5169992) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1285308) q[1];
sx q[1];
rz(-1.3088041) q[1];
sx q[1];
rz(-0.2688516) q[1];
rz(-0.92658073) q[3];
sx q[3];
rz(-0.80873064) q[3];
sx q[3];
rz(-0.96264983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.076685585) q[2];
sx q[2];
rz(-1.6904597) q[2];
sx q[2];
rz(-1.1265075) q[2];
rz(1.7846151) q[3];
sx q[3];
rz(-0.49961909) q[3];
sx q[3];
rz(-0.13993851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5821238) q[0];
sx q[0];
rz(-2.421565) q[0];
sx q[0];
rz(-0.35991392) q[0];
rz(0.24078807) q[1];
sx q[1];
rz(-1.9937932) q[1];
sx q[1];
rz(-0.37318939) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8561962) q[0];
sx q[0];
rz(-1.4715949) q[0];
sx q[0];
rz(-1.7147934) q[0];
rz(-pi) q[1];
rz(1.3360489) q[2];
sx q[2];
rz(-2.1214607) q[2];
sx q[2];
rz(2.5900813) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2265985) q[1];
sx q[1];
rz(-1.0553331) q[1];
sx q[1];
rz(1.9644009) q[1];
rz(-2.7591428) q[3];
sx q[3];
rz(-2.6646864) q[3];
sx q[3];
rz(-2.675569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.084426247) q[2];
sx q[2];
rz(-1.7095704) q[2];
sx q[2];
rz(2.2331494) q[2];
rz(-1.069979) q[3];
sx q[3];
rz(-1.1845651) q[3];
sx q[3];
rz(2.158304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.95803607) q[0];
sx q[0];
rz(-0.89360845) q[0];
sx q[0];
rz(0.91304427) q[0];
rz(-0.68963447) q[1];
sx q[1];
rz(-0.90877405) q[1];
sx q[1];
rz(-2.3675809) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.000074) q[0];
sx q[0];
rz(-0.98709269) q[0];
sx q[0];
rz(1.0144177) q[0];
rz(-pi) q[1];
x q[1];
rz(0.795587) q[2];
sx q[2];
rz(-1.5485816) q[2];
sx q[2];
rz(-1.8088248) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4178631) q[1];
sx q[1];
rz(-2.6681719) q[1];
sx q[1];
rz(1.141505) q[1];
rz(1.6004531) q[3];
sx q[3];
rz(-1.0231442) q[3];
sx q[3];
rz(0.95049196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.245605) q[2];
sx q[2];
rz(-1.3934803) q[2];
sx q[2];
rz(0.73949933) q[2];
rz(-1.5064404) q[3];
sx q[3];
rz(-2.2061901) q[3];
sx q[3];
rz(-0.80772775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6226115) q[0];
sx q[0];
rz(-2.3464572) q[0];
sx q[0];
rz(-1.125289) q[0];
rz(-3.0033424) q[1];
sx q[1];
rz(-1.6743276) q[1];
sx q[1];
rz(-1.5672055) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0042717) q[0];
sx q[0];
rz(-1.204365) q[0];
sx q[0];
rz(1.9204813) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2786516) q[2];
sx q[2];
rz(-1.9685479) q[2];
sx q[2];
rz(-1.5270555) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.36260819) q[1];
sx q[1];
rz(-1.6951188) q[1];
sx q[1];
rz(-0.26047867) q[1];
rz(2.8254824) q[3];
sx q[3];
rz(-0.5721285) q[3];
sx q[3];
rz(0.24219777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.77341998) q[2];
sx q[2];
rz(-2.07351) q[2];
sx q[2];
rz(2.0424021) q[2];
rz(0.98569551) q[3];
sx q[3];
rz(-1.6253977) q[3];
sx q[3];
rz(-2.709008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2987591) q[0];
sx q[0];
rz(-2.1465813) q[0];
sx q[0];
rz(-2.3327935) q[0];
rz(2.47593) q[1];
sx q[1];
rz(-0.97949615) q[1];
sx q[1];
rz(2.1601802) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4047667) q[0];
sx q[0];
rz(-0.5739218) q[0];
sx q[0];
rz(1.4653652) q[0];
rz(-pi) q[1];
rz(0.50534525) q[2];
sx q[2];
rz(-1.7957558) q[2];
sx q[2];
rz(-2.7749429) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8725207) q[1];
sx q[1];
rz(-1.4351294) q[1];
sx q[1];
rz(-1.1024464) q[1];
rz(-pi) q[2];
rz(3.0043169) q[3];
sx q[3];
rz(-2.2386754) q[3];
sx q[3];
rz(-1.5811416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.064528) q[2];
sx q[2];
rz(-1.6669824) q[2];
sx q[2];
rz(-2.7725753) q[2];
rz(-1.4024233) q[3];
sx q[3];
rz(-2.2544421) q[3];
sx q[3];
rz(1.3212475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8333261) q[0];
sx q[0];
rz(-2.3127191) q[0];
sx q[0];
rz(-2.8114124) q[0];
rz(0.45686832) q[1];
sx q[1];
rz(-2.7062682) q[1];
sx q[1];
rz(-0.59828573) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9015181) q[0];
sx q[0];
rz(-0.80722729) q[0];
sx q[0];
rz(-0.62348311) q[0];
rz(1.4315375) q[2];
sx q[2];
rz(-2.6684847) q[2];
sx q[2];
rz(2.3654761) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.51914267) q[1];
sx q[1];
rz(-2.7323117) q[1];
sx q[1];
rz(2.5342788) q[1];
x q[2];
rz(0.98294799) q[3];
sx q[3];
rz(-1.0025257) q[3];
sx q[3];
rz(-2.1592888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2716219) q[2];
sx q[2];
rz(-2.2651256) q[2];
sx q[2];
rz(-0.064621933) q[2];
rz(-1.6914852) q[3];
sx q[3];
rz(-2.7414069) q[3];
sx q[3];
rz(1.5959285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86218631) q[0];
sx q[0];
rz(-0.33410826) q[0];
sx q[0];
rz(-2.620328) q[0];
rz(1.0598671) q[1];
sx q[1];
rz(-1.1618549) q[1];
sx q[1];
rz(2.1102139) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8169075) q[0];
sx q[0];
rz(-0.74231901) q[0];
sx q[0];
rz(-3.0317505) q[0];
x q[1];
rz(0.29287405) q[2];
sx q[2];
rz(-1.4559064) q[2];
sx q[2];
rz(3.1174768) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0963991) q[1];
sx q[1];
rz(-0.61122433) q[1];
sx q[1];
rz(1.0525714) q[1];
rz(-pi) q[2];
rz(-1.6691339) q[3];
sx q[3];
rz(-1.0957484) q[3];
sx q[3];
rz(2.6779384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3970268) q[2];
sx q[2];
rz(-0.72781813) q[2];
sx q[2];
rz(-3.0542206) q[2];
rz(1.9258063) q[3];
sx q[3];
rz(-1.4601424) q[3];
sx q[3];
rz(2.9714835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66607296) q[0];
sx q[0];
rz(-2.2505794) q[0];
sx q[0];
rz(-1.1391621) q[0];
rz(-2.6423404) q[1];
sx q[1];
rz(-1.4491932) q[1];
sx q[1];
rz(1.8081236) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77908726) q[0];
sx q[0];
rz(-1.6769857) q[0];
sx q[0];
rz(-1.6615191) q[0];
x q[1];
rz(3.0917909) q[2];
sx q[2];
rz(-1.8254455) q[2];
sx q[2];
rz(-1.7083502) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1386316) q[1];
sx q[1];
rz(-1.8233144) q[1];
sx q[1];
rz(1.3481524) q[1];
x q[2];
rz(-0.91657775) q[3];
sx q[3];
rz(-0.98482705) q[3];
sx q[3];
rz(0.41205349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6405876) q[2];
sx q[2];
rz(-1.2747526) q[2];
sx q[2];
rz(0.28323832) q[2];
rz(-1.9404274) q[3];
sx q[3];
rz(-0.675942) q[3];
sx q[3];
rz(2.0961608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82861154) q[0];
sx q[0];
rz(-2.3647478) q[0];
sx q[0];
rz(1.1051529) q[0];
rz(3.037187) q[1];
sx q[1];
rz(-1.5515635) q[1];
sx q[1];
rz(-1.9952231) q[1];
rz(2.636006) q[2];
sx q[2];
rz(-1.7442295) q[2];
sx q[2];
rz(-2.8019047) q[2];
rz(0.66815175) q[3];
sx q[3];
rz(-1.2378319) q[3];
sx q[3];
rz(1.1938865) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
