OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0661434) q[0];
sx q[0];
rz(-2.0976522) q[0];
sx q[0];
rz(3.1312842) q[0];
rz(2.161624) q[1];
sx q[1];
rz(-1.4449395) q[1];
sx q[1];
rz(2.565032) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.071237) q[0];
sx q[0];
rz(-2.8094387) q[0];
sx q[0];
rz(2.7624056) q[0];
x q[1];
rz(-2.2498807) q[2];
sx q[2];
rz(-1.0701961) q[2];
sx q[2];
rz(-2.5426189) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.7652055) q[1];
sx q[1];
rz(-2.58316) q[1];
sx q[1];
rz(-2.8064578) q[1];
x q[2];
rz(2.303894) q[3];
sx q[3];
rz(-1.4464966) q[3];
sx q[3];
rz(-1.705114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.52446857) q[2];
sx q[2];
rz(-1.4602129) q[2];
sx q[2];
rz(-1.2855533) q[2];
rz(1.5517976) q[3];
sx q[3];
rz(-1.0714622) q[3];
sx q[3];
rz(-2.0901399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1262421) q[0];
sx q[0];
rz(-1.9168357) q[0];
sx q[0];
rz(-2.8958564) q[0];
rz(-1.0579146) q[1];
sx q[1];
rz(-1.825288) q[1];
sx q[1];
rz(2.8071075) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4870105) q[0];
sx q[0];
rz(-2.062487) q[0];
sx q[0];
rz(-2.0915178) q[0];
x q[1];
rz(-1.9742786) q[2];
sx q[2];
rz(-0.85131391) q[2];
sx q[2];
rz(3.1297504) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.04131728) q[1];
sx q[1];
rz(-1.9875437) q[1];
sx q[1];
rz(-2.7303425) q[1];
rz(0.40615079) q[3];
sx q[3];
rz(-0.63780071) q[3];
sx q[3];
rz(-2.8075308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1841396) q[2];
sx q[2];
rz(-0.89749557) q[2];
sx q[2];
rz(-1.3857566) q[2];
rz(0.98006788) q[3];
sx q[3];
rz(-2.6762784) q[3];
sx q[3];
rz(-0.050962713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7053213) q[0];
sx q[0];
rz(-2.1846117) q[0];
sx q[0];
rz(1.7514239) q[0];
rz(2.7888489) q[1];
sx q[1];
rz(-0.91068641) q[1];
sx q[1];
rz(-2.5211451) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13201319) q[0];
sx q[0];
rz(-1.4963829) q[0];
sx q[0];
rz(-0.090601765) q[0];
rz(3.0376126) q[2];
sx q[2];
rz(-2.1562088) q[2];
sx q[2];
rz(-0.81987655) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7812278) q[1];
sx q[1];
rz(-1.8958127) q[1];
sx q[1];
rz(-2.8311391) q[1];
x q[2];
rz(-2.6246715) q[3];
sx q[3];
rz(-0.94603387) q[3];
sx q[3];
rz(2.3464835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0188401) q[2];
sx q[2];
rz(-2.7285125) q[2];
sx q[2];
rz(1.2472461) q[2];
rz(-2.9774169) q[3];
sx q[3];
rz(-1.434606) q[3];
sx q[3];
rz(0.629614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71488798) q[0];
sx q[0];
rz(-2.1737104) q[0];
sx q[0];
rz(-1.2870652) q[0];
rz(2.5413068) q[1];
sx q[1];
rz(-1.3515819) q[1];
sx q[1];
rz(0.90369019) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1290734) q[0];
sx q[0];
rz(-3.0154069) q[0];
sx q[0];
rz(0.55380765) q[0];
rz(-pi) q[1];
rz(2.972347) q[2];
sx q[2];
rz(-2.7566559) q[2];
sx q[2];
rz(0.16916179) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7052482) q[1];
sx q[1];
rz(-1.7503993) q[1];
sx q[1];
rz(-0.257538) q[1];
x q[2];
rz(-1.8925531) q[3];
sx q[3];
rz(-1.9638954) q[3];
sx q[3];
rz(-2.2571795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6220182) q[2];
sx q[2];
rz(-0.833424) q[2];
sx q[2];
rz(-0.77345094) q[2];
rz(1.8526239) q[3];
sx q[3];
rz(-0.52923146) q[3];
sx q[3];
rz(2.4066063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89161038) q[0];
sx q[0];
rz(-2.1365428) q[0];
sx q[0];
rz(-2.6494359) q[0];
rz(-2.6026169) q[1];
sx q[1];
rz(-1.4424126) q[1];
sx q[1];
rz(-1.0999365) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.625222) q[0];
sx q[0];
rz(-1.3110135) q[0];
sx q[0];
rz(2.7894839) q[0];
x q[1];
rz(-0.01077588) q[2];
sx q[2];
rz(-2.1824565) q[2];
sx q[2];
rz(1.3792559) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.089842794) q[1];
sx q[1];
rz(-0.81498346) q[1];
sx q[1];
rz(-3.0137193) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8995011) q[3];
sx q[3];
rz(-1.1100169) q[3];
sx q[3];
rz(-2.8976909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.34053549) q[2];
sx q[2];
rz(-0.33581442) q[2];
sx q[2];
rz(-0.56279969) q[2];
rz(0.69989145) q[3];
sx q[3];
rz(-1.2908582) q[3];
sx q[3];
rz(-2.8904397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7971147) q[0];
sx q[0];
rz(-0.7239224) q[0];
sx q[0];
rz(-2.7864454) q[0];
rz(-1.9225325) q[1];
sx q[1];
rz(-1.2390169) q[1];
sx q[1];
rz(-2.6893137) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1817976) q[0];
sx q[0];
rz(-0.99645146) q[0];
sx q[0];
rz(2.6622165) q[0];
x q[1];
rz(-2.7554871) q[2];
sx q[2];
rz(-1.9506644) q[2];
sx q[2];
rz(-2.7136193) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.86231316) q[1];
sx q[1];
rz(-1.8932027) q[1];
sx q[1];
rz(-0.025397852) q[1];
x q[2];
rz(-2.5127453) q[3];
sx q[3];
rz(-1.4553242) q[3];
sx q[3];
rz(-0.94146282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.63048116) q[2];
sx q[2];
rz(-0.52383542) q[2];
sx q[2];
rz(0.13709489) q[2];
rz(1.5824205) q[3];
sx q[3];
rz(-1.987792) q[3];
sx q[3];
rz(-2.3649575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(0.95098507) q[0];
sx q[0];
rz(-0.76222104) q[0];
sx q[0];
rz(1.5964339) q[0];
rz(0.51721382) q[1];
sx q[1];
rz(-1.1511753) q[1];
sx q[1];
rz(-2.1545765) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2632946) q[0];
sx q[0];
rz(-1.7016101) q[0];
sx q[0];
rz(1.592171) q[0];
rz(-pi) q[1];
rz(0.9632684) q[2];
sx q[2];
rz(-1.2973366) q[2];
sx q[2];
rz(-0.69404049) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2199041) q[1];
sx q[1];
rz(-2.7877586) q[1];
sx q[1];
rz(-0.71531957) q[1];
rz(-pi) q[2];
rz(-1.889398) q[3];
sx q[3];
rz(-1.4907537) q[3];
sx q[3];
rz(-2.6032676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.78642693) q[2];
sx q[2];
rz(-1.0309018) q[2];
sx q[2];
rz(0.008358566) q[2];
rz(0.23056325) q[3];
sx q[3];
rz(-1.9356666) q[3];
sx q[3];
rz(-1.4768538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.01345988) q[0];
sx q[0];
rz(-3.1210493) q[0];
sx q[0];
rz(-1.531456) q[0];
rz(-1.6186591) q[1];
sx q[1];
rz(-1.5956655) q[1];
sx q[1];
rz(-1.5628372) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7195936) q[0];
sx q[0];
rz(-2.5223456) q[0];
sx q[0];
rz(-2.8887755) q[0];
rz(-pi) q[1];
rz(-0.35308102) q[2];
sx q[2];
rz(-0.39829474) q[2];
sx q[2];
rz(-2.9495267) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2819689) q[1];
sx q[1];
rz(-2.0323557) q[1];
sx q[1];
rz(2.0412316) q[1];
x q[2];
rz(-1.8914901) q[3];
sx q[3];
rz(-0.9281635) q[3];
sx q[3];
rz(-1.927141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.83850399) q[2];
sx q[2];
rz(-1.9946626) q[2];
sx q[2];
rz(0.037503555) q[2];
rz(2.904902) q[3];
sx q[3];
rz(-2.762837) q[3];
sx q[3];
rz(-1.8481351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8732052) q[0];
sx q[0];
rz(-2.3516042) q[0];
sx q[0];
rz(1.0733806) q[0];
rz(2.7997596) q[1];
sx q[1];
rz(-2.5624202) q[1];
sx q[1];
rz(2.5198708) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8499924) q[0];
sx q[0];
rz(-1.9555641) q[0];
sx q[0];
rz(-1.8509393) q[0];
x q[1];
rz(0.15839496) q[2];
sx q[2];
rz(-2.4753053) q[2];
sx q[2];
rz(-1.0384384) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7555862) q[1];
sx q[1];
rz(-1.2666128) q[1];
sx q[1];
rz(-2.6261397) q[1];
x q[2];
rz(3.000876) q[3];
sx q[3];
rz(-1.4840115) q[3];
sx q[3];
rz(-1.113689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6217893) q[2];
sx q[2];
rz(-0.18814627) q[2];
sx q[2];
rz(-0.56358799) q[2];
rz(-1.8300736) q[3];
sx q[3];
rz(-1.2389641) q[3];
sx q[3];
rz(-0.81609503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(2.9594864) q[0];
sx q[0];
rz(-0.86903787) q[0];
sx q[0];
rz(0.8748138) q[0];
rz(1.733571) q[1];
sx q[1];
rz(-0.66649109) q[1];
sx q[1];
rz(0.63546884) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4548774) q[0];
sx q[0];
rz(-1.8679163) q[0];
sx q[0];
rz(2.3700506) q[0];
x q[1];
rz(-2.1383912) q[2];
sx q[2];
rz(-0.48241189) q[2];
sx q[2];
rz(0.29200992) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8567896) q[1];
sx q[1];
rz(-1.9724047) q[1];
sx q[1];
rz(2.1339586) q[1];
rz(-pi) q[2];
rz(-2.3395679) q[3];
sx q[3];
rz(-2.4200508) q[3];
sx q[3];
rz(-2.3156543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1982939) q[2];
sx q[2];
rz(-2.0512927) q[2];
sx q[2];
rz(-2.3401006) q[2];
rz(-1.9258026) q[3];
sx q[3];
rz(-0.91167584) q[3];
sx q[3];
rz(0.041778684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53032482) q[0];
sx q[0];
rz(-1.4401191) q[0];
sx q[0];
rz(-2.8339207) q[0];
rz(1.8511741) q[1];
sx q[1];
rz(-2.1967874) q[1];
sx q[1];
rz(-2.8312942) q[1];
rz(0.18193131) q[2];
sx q[2];
rz(-1.7040737) q[2];
sx q[2];
rz(-2.2752442) q[2];
rz(1.7054059) q[3];
sx q[3];
rz(-2.2231839) q[3];
sx q[3];
rz(-0.25683944) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
