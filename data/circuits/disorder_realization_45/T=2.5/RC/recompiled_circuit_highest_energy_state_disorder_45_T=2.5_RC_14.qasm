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
rz(-0.51195872) q[0];
sx q[0];
rz(-2.7380138) q[0];
sx q[0];
rz(3.1374748) q[0];
rz(0.68139684) q[1];
sx q[1];
rz(-1.2702785) q[1];
sx q[1];
rz(-1.7008002) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3128634) q[0];
sx q[0];
rz(-1.848319) q[0];
sx q[0];
rz(2.9636431) q[0];
rz(1.2295159) q[2];
sx q[2];
rz(-0.84898872) q[2];
sx q[2];
rz(-2.316458) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5087252) q[1];
sx q[1];
rz(-2.3687216) q[1];
sx q[1];
rz(0.51846026) q[1];
rz(-pi) q[2];
rz(-0.012083455) q[3];
sx q[3];
rz(-1.3971431) q[3];
sx q[3];
rz(2.1326667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9422841) q[2];
sx q[2];
rz(-2.0950623) q[2];
sx q[2];
rz(-2.671177) q[2];
rz(2.2453902) q[3];
sx q[3];
rz(-1.4678518) q[3];
sx q[3];
rz(1.5907653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6386221) q[0];
sx q[0];
rz(-2.6574385) q[0];
sx q[0];
rz(0.36889398) q[0];
rz(-2.6781354) q[1];
sx q[1];
rz(-1.2838485) q[1];
sx q[1];
rz(2.8772112) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3437816) q[0];
sx q[0];
rz(-2.4759295) q[0];
sx q[0];
rz(1.8209609) q[0];
rz(-pi) q[1];
rz(0.55230852) q[2];
sx q[2];
rz(-2.8495516) q[2];
sx q[2];
rz(0.24667443) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0708587) q[1];
sx q[1];
rz(-1.8740591) q[1];
sx q[1];
rz(1.2649078) q[1];
x q[2];
rz(-0.12118487) q[3];
sx q[3];
rz(-0.78599343) q[3];
sx q[3];
rz(-1.3898897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.6139063) q[2];
sx q[2];
rz(-0.1656342) q[2];
sx q[2];
rz(-1.3045093) q[2];
rz(-0.3197318) q[3];
sx q[3];
rz(-0.78800646) q[3];
sx q[3];
rz(0.50362292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17427915) q[0];
sx q[0];
rz(-0.19198424) q[0];
sx q[0];
rz(2.018003) q[0];
rz(-2.7536821) q[1];
sx q[1];
rz(-1.265641) q[1];
sx q[1];
rz(0.33043114) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.743227) q[0];
sx q[0];
rz(-1.3083698) q[0];
sx q[0];
rz(-2.9578379) q[0];
x q[1];
rz(-0.60744384) q[2];
sx q[2];
rz(-1.6472378) q[2];
sx q[2];
rz(-2.3747403) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.36281113) q[1];
sx q[1];
rz(-1.2515725) q[1];
sx q[1];
rz(-0.37324639) q[1];
x q[2];
rz(-0.64066842) q[3];
sx q[3];
rz(-2.0694149) q[3];
sx q[3];
rz(-0.22646587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.72184163) q[2];
sx q[2];
rz(-1.9998877) q[2];
sx q[2];
rz(1.9640131) q[2];
rz(-1.8925331) q[3];
sx q[3];
rz(-1.8355337) q[3];
sx q[3];
rz(-0.038399847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.053442001) q[0];
sx q[0];
rz(-1.4295239) q[0];
sx q[0];
rz(0.088726774) q[0];
rz(0.42452043) q[1];
sx q[1];
rz(-2.4219234) q[1];
sx q[1];
rz(-1.7321865) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8363894) q[0];
sx q[0];
rz(-0.38581784) q[0];
sx q[0];
rz(1.057339) q[0];
rz(-pi) q[1];
rz(2.2327044) q[2];
sx q[2];
rz(-0.7134921) q[2];
sx q[2];
rz(-2.5475218) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.030050412) q[1];
sx q[1];
rz(-0.71978986) q[1];
sx q[1];
rz(-0.58776249) q[1];
rz(-pi) q[2];
rz(-0.32281422) q[3];
sx q[3];
rz(-1.0720616) q[3];
sx q[3];
rz(2.7018967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.38922629) q[2];
sx q[2];
rz(-1.0800635) q[2];
sx q[2];
rz(-2.7092095) q[2];
rz(-2.6835175) q[3];
sx q[3];
rz(-1.4832486) q[3];
sx q[3];
rz(2.1336011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78419375) q[0];
sx q[0];
rz(-2.9965239) q[0];
sx q[0];
rz(2.8302637) q[0];
rz(-1.1773102) q[1];
sx q[1];
rz(-1.9639587) q[1];
sx q[1];
rz(2.6572773) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70894079) q[0];
sx q[0];
rz(-1.4756894) q[0];
sx q[0];
rz(1.7300138) q[0];
rz(-pi) q[1];
rz(1.2782477) q[2];
sx q[2];
rz(-1.6451837) q[2];
sx q[2];
rz(2.0015581) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2023485) q[1];
sx q[1];
rz(-2.7396297) q[1];
sx q[1];
rz(-1.9415929) q[1];
rz(-pi) q[2];
rz(2.9919162) q[3];
sx q[3];
rz(-1.6972741) q[3];
sx q[3];
rz(1.5095131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3207265) q[2];
sx q[2];
rz(-1.4297239) q[2];
sx q[2];
rz(2.7166264) q[2];
rz(3.037437) q[3];
sx q[3];
rz(-2.7724373) q[3];
sx q[3];
rz(0.66515508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3686309) q[0];
sx q[0];
rz(-0.63987982) q[0];
sx q[0];
rz(-2.9848918) q[0];
rz(-0.26198298) q[1];
sx q[1];
rz(-1.9888839) q[1];
sx q[1];
rz(1.6798457) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7279434) q[0];
sx q[0];
rz(-0.62725146) q[0];
sx q[0];
rz(1.2722737) q[0];
rz(1.842672) q[2];
sx q[2];
rz(-0.71496925) q[2];
sx q[2];
rz(0.78503099) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.78292978) q[1];
sx q[1];
rz(-1.6921669) q[1];
sx q[1];
rz(0.18106457) q[1];
rz(-pi) q[2];
rz(-2.7191592) q[3];
sx q[3];
rz(-2.3452873) q[3];
sx q[3];
rz(0.56047201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.67019749) q[2];
sx q[2];
rz(-1.6452226) q[2];
sx q[2];
rz(0.71693286) q[2];
rz(0.21643058) q[3];
sx q[3];
rz(-1.2856893) q[3];
sx q[3];
rz(-1.1855679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4099429) q[0];
sx q[0];
rz(-0.2922295) q[0];
sx q[0];
rz(1.8753847) q[0];
rz(0.75824291) q[1];
sx q[1];
rz(-1.1781324) q[1];
sx q[1];
rz(-2.0936802) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78572014) q[0];
sx q[0];
rz(-1.2730755) q[0];
sx q[0];
rz(-2.1773318) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1282461) q[2];
sx q[2];
rz(-1.1801022) q[2];
sx q[2];
rz(1.3732131) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0430773) q[1];
sx q[1];
rz(-1.2081573) q[1];
sx q[1];
rz(0.074438268) q[1];
x q[2];
rz(1.8147477) q[3];
sx q[3];
rz(-2.0097964) q[3];
sx q[3];
rz(1.148996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8083814) q[2];
sx q[2];
rz(-2.1344678) q[2];
sx q[2];
rz(-3.132931) q[2];
rz(-1.0041142) q[3];
sx q[3];
rz(-2.6486371) q[3];
sx q[3];
rz(-2.139411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36326161) q[0];
sx q[0];
rz(-0.14334981) q[0];
sx q[0];
rz(-0.95247254) q[0];
rz(1.5337503) q[1];
sx q[1];
rz(-1.8550823) q[1];
sx q[1];
rz(1.0362157) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0909611) q[0];
sx q[0];
rz(-1.0337752) q[0];
sx q[0];
rz(-1.0616674) q[0];
x q[1];
rz(2.6626946) q[2];
sx q[2];
rz(-0.69456813) q[2];
sx q[2];
rz(-2.640129) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.10401512) q[1];
sx q[1];
rz(-1.0733782) q[1];
sx q[1];
rz(-0.90105547) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0228859) q[3];
sx q[3];
rz(-2.2452659) q[3];
sx q[3];
rz(-1.7728286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8448392) q[2];
sx q[2];
rz(-0.71594683) q[2];
sx q[2];
rz(2.1412264) q[2];
rz(-1.865546) q[3];
sx q[3];
rz(-1.4714656) q[3];
sx q[3];
rz(-2.0744417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.100383) q[0];
sx q[0];
rz(-0.90827933) q[0];
sx q[0];
rz(0.25417438) q[0];
rz(1.8036448) q[1];
sx q[1];
rz(-0.86054069) q[1];
sx q[1];
rz(2.3635704) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7800919) q[0];
sx q[0];
rz(-2.2169212) q[0];
sx q[0];
rz(1.161399) q[0];
rz(2.940321) q[2];
sx q[2];
rz(-1.4106371) q[2];
sx q[2];
rz(1.8294272) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.48315037) q[1];
sx q[1];
rz(-1.1643693) q[1];
sx q[1];
rz(-2.5220925) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.098962231) q[3];
sx q[3];
rz(-1.785019) q[3];
sx q[3];
rz(0.19161397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0622056) q[2];
sx q[2];
rz(-1.5176682) q[2];
sx q[2];
rz(0.13266955) q[2];
rz(2.0343871) q[3];
sx q[3];
rz(-2.5532494) q[3];
sx q[3];
rz(1.4440822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74141136) q[0];
sx q[0];
rz(-1.1962698) q[0];
sx q[0];
rz(-0.74914002) q[0];
rz(-2.6465042) q[1];
sx q[1];
rz(-1.2352713) q[1];
sx q[1];
rz(-1.3912158) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.036450841) q[0];
sx q[0];
rz(-0.81691475) q[0];
sx q[0];
rz(-0.63936279) q[0];
x q[1];
rz(-0.40410186) q[2];
sx q[2];
rz(-1.1617336) q[2];
sx q[2];
rz(-2.7210195) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.039469214) q[1];
sx q[1];
rz(-1.9774924) q[1];
sx q[1];
rz(1.6182093) q[1];
rz(-2.0999317) q[3];
sx q[3];
rz(-1.1253998) q[3];
sx q[3];
rz(-0.22280773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.64080015) q[2];
sx q[2];
rz(-1.6668789) q[2];
sx q[2];
rz(2.2643209) q[2];
rz(2.6194465) q[3];
sx q[3];
rz(-0.79972655) q[3];
sx q[3];
rz(-1.6845901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40182879) q[0];
sx q[0];
rz(-0.54720989) q[0];
sx q[0];
rz(-1.3450958) q[0];
rz(1.8632035) q[1];
sx q[1];
rz(-0.31247333) q[1];
sx q[1];
rz(-0.023275274) q[1];
rz(-1.5398434) q[2];
sx q[2];
rz(-2.2014115) q[2];
sx q[2];
rz(1.9104107) q[2];
rz(0.69576453) q[3];
sx q[3];
rz(-2.5740665) q[3];
sx q[3];
rz(3.1181581) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
