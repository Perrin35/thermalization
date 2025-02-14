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
rz(-1.902154) q[0];
sx q[0];
rz(-1.3286989) q[0];
sx q[0];
rz(-2.9297096) q[0];
rz(1.7243241) q[1];
sx q[1];
rz(-0.53242004) q[1];
sx q[1];
rz(-0.37791696) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.657562) q[0];
sx q[0];
rz(-3.1320509) q[0];
sx q[0];
rz(-1.6173167) q[0];
x q[1];
rz(1.9340408) q[2];
sx q[2];
rz(-0.9245199) q[2];
sx q[2];
rz(0.56528795) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.35740543) q[1];
sx q[1];
rz(-0.94036803) q[1];
sx q[1];
rz(-1.9131843) q[1];
rz(-2.3842281) q[3];
sx q[3];
rz(-0.93776449) q[3];
sx q[3];
rz(-0.6642793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2892896) q[2];
sx q[2];
rz(-1.0785582) q[2];
sx q[2];
rz(-0.079785384) q[2];
rz(-0.96528178) q[3];
sx q[3];
rz(-1.4502757) q[3];
sx q[3];
rz(0.7307581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47736436) q[0];
sx q[0];
rz(-1.0034765) q[0];
sx q[0];
rz(3.0526414) q[0];
rz(1.3720007) q[1];
sx q[1];
rz(-1.4274495) q[1];
sx q[1];
rz(-1.1044097) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7331284) q[0];
sx q[0];
rz(-0.95703546) q[0];
sx q[0];
rz(-1.2888786) q[0];
rz(-1.0816108) q[2];
sx q[2];
rz(-2.2913165) q[2];
sx q[2];
rz(-0.94948506) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4648393) q[1];
sx q[1];
rz(-2.0130806) q[1];
sx q[1];
rz(2.9763336) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6703963) q[3];
sx q[3];
rz(-2.1079113) q[3];
sx q[3];
rz(-0.84505862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8770807) q[2];
sx q[2];
rz(-1.1254213) q[2];
sx q[2];
rz(1.4667) q[2];
rz(-3.1125715) q[3];
sx q[3];
rz(-1.0163739) q[3];
sx q[3];
rz(-0.69028729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90917176) q[0];
sx q[0];
rz(-0.066815289) q[0];
sx q[0];
rz(0.27798852) q[0];
rz(-1.4311283) q[1];
sx q[1];
rz(-0.93228308) q[1];
sx q[1];
rz(-2.7412282) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6872266) q[0];
sx q[0];
rz(-2.3355744) q[0];
sx q[0];
rz(0.63712742) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.68636559) q[2];
sx q[2];
rz(-1.371742) q[2];
sx q[2];
rz(0.65716568) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4440205) q[1];
sx q[1];
rz(-2.6409147) q[1];
sx q[1];
rz(2.8050007) q[1];
rz(-2.8002162) q[3];
sx q[3];
rz(-1.8601283) q[3];
sx q[3];
rz(-0.067726243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6843188) q[2];
sx q[2];
rz(-1.6088586) q[2];
sx q[2];
rz(-2.686783) q[2];
rz(-1.8999892) q[3];
sx q[3];
rz(-0.95507115) q[3];
sx q[3];
rz(0.72030592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1860745) q[0];
sx q[0];
rz(-1.8121413) q[0];
sx q[0];
rz(3.1410826) q[0];
rz(-0.60091248) q[1];
sx q[1];
rz(-0.83952236) q[1];
sx q[1];
rz(-2.9972163) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6116981) q[0];
sx q[0];
rz(-1.3641883) q[0];
sx q[0];
rz(1.0103465) q[0];
rz(-pi) q[1];
x q[1];
rz(0.52835502) q[2];
sx q[2];
rz(-0.29137416) q[2];
sx q[2];
rz(-1.0519415) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2541009) q[1];
sx q[1];
rz(-2.2057167) q[1];
sx q[1];
rz(-1.8545521) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9668429) q[3];
sx q[3];
rz(-0.55395444) q[3];
sx q[3];
rz(2.9308211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.942975) q[2];
sx q[2];
rz(-1.4451507) q[2];
sx q[2];
rz(0.96251881) q[2];
rz(1.6019542) q[3];
sx q[3];
rz(-1.736015) q[3];
sx q[3];
rz(-0.34077728) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2365504) q[0];
sx q[0];
rz(-1.4555229) q[0];
sx q[0];
rz(-1.9566253) q[0];
rz(0.22625893) q[1];
sx q[1];
rz(-0.87892756) q[1];
sx q[1];
rz(1.0386946) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13030355) q[0];
sx q[0];
rz(-1.652392) q[0];
sx q[0];
rz(2.5528583) q[0];
rz(2.4034385) q[2];
sx q[2];
rz(-2.0655491) q[2];
sx q[2];
rz(0.02502266) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9925527) q[1];
sx q[1];
rz(-2.2254308) q[1];
sx q[1];
rz(2.4813985) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0108972) q[3];
sx q[3];
rz(-1.6301042) q[3];
sx q[3];
rz(1.9594994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.431939) q[2];
sx q[2];
rz(-2.4434872) q[2];
sx q[2];
rz(-0.31614885) q[2];
rz(1.4656674) q[3];
sx q[3];
rz(-1.1301872) q[3];
sx q[3];
rz(-1.7620311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50452152) q[0];
sx q[0];
rz(-2.2383454) q[0];
sx q[0];
rz(-2.0380518) q[0];
rz(0.48031131) q[1];
sx q[1];
rz(-0.66245285) q[1];
sx q[1];
rz(-2.5441817) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6107268) q[0];
sx q[0];
rz(-1.382575) q[0];
sx q[0];
rz(-1.782531) q[0];
rz(-pi) q[1];
rz(-0.2603724) q[2];
sx q[2];
rz(-0.7625167) q[2];
sx q[2];
rz(2.0338361) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.5525517) q[1];
sx q[1];
rz(-1.8218166) q[1];
sx q[1];
rz(-1.4000721) q[1];
x q[2];
rz(1.1998981) q[3];
sx q[3];
rz(-1.6565595) q[3];
sx q[3];
rz(-0.24468064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9395113) q[2];
sx q[2];
rz(-1.5152405) q[2];
sx q[2];
rz(-1.0763947) q[2];
rz(1.1427897) q[3];
sx q[3];
rz(-2.3484774) q[3];
sx q[3];
rz(0.49016652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.48653212) q[0];
sx q[0];
rz(-1.158411) q[0];
sx q[0];
rz(3.1031188) q[0];
rz(3.0768652) q[1];
sx q[1];
rz(-1.7579301) q[1];
sx q[1];
rz(2.9077392) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42002871) q[0];
sx q[0];
rz(-1.7917243) q[0];
sx q[0];
rz(1.7608282) q[0];
rz(-2.8296956) q[2];
sx q[2];
rz(-1.9210749) q[2];
sx q[2];
rz(0.52679449) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5906787) q[1];
sx q[1];
rz(-1.4812902) q[1];
sx q[1];
rz(-2.4924335) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7247808) q[3];
sx q[3];
rz(-0.14188611) q[3];
sx q[3];
rz(-2.8519423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6528299) q[2];
sx q[2];
rz(-1.5583928) q[2];
sx q[2];
rz(2.5433507) q[2];
rz(0.11387842) q[3];
sx q[3];
rz(-1.3860393) q[3];
sx q[3];
rz(-0.84754506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4050196) q[0];
sx q[0];
rz(-2.2990655) q[0];
sx q[0];
rz(2.8952428) q[0];
rz(-1.8440638) q[1];
sx q[1];
rz(-2.0228701) q[1];
sx q[1];
rz(-2.6447703) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9333736) q[0];
sx q[0];
rz(-0.95666203) q[0];
sx q[0];
rz(1.0557515) q[0];
rz(-pi) q[1];
x q[1];
rz(0.17077568) q[2];
sx q[2];
rz(-0.71238067) q[2];
sx q[2];
rz(0.26584372) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5766746) q[1];
sx q[1];
rz(-1.0905301) q[1];
sx q[1];
rz(0.26580055) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0629749) q[3];
sx q[3];
rz(-1.7075305) q[3];
sx q[3];
rz(1.1059974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7184427) q[2];
sx q[2];
rz(-0.53360525) q[2];
sx q[2];
rz(1.144484) q[2];
rz(1.1540958) q[3];
sx q[3];
rz(-1.7172979) q[3];
sx q[3];
rz(-0.19788876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6941187) q[0];
sx q[0];
rz(-0.23451528) q[0];
sx q[0];
rz(3.0754454) q[0];
rz(-1.1478395) q[1];
sx q[1];
rz(-1.3775974) q[1];
sx q[1];
rz(0.60417169) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9420931) q[0];
sx q[0];
rz(-1.2954933) q[0];
sx q[0];
rz(-1.784174) q[0];
x q[1];
rz(-0.2653052) q[2];
sx q[2];
rz(-0.54232208) q[2];
sx q[2];
rz(0.10077439) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7614038) q[1];
sx q[1];
rz(-1.0127002) q[1];
sx q[1];
rz(1.3437043) q[1];
rz(-2.3202032) q[3];
sx q[3];
rz(-0.8474955) q[3];
sx q[3];
rz(2.5625474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8784647) q[2];
sx q[2];
rz(-0.85004127) q[2];
sx q[2];
rz(-2.7086332) q[2];
rz(-1.2290907) q[3];
sx q[3];
rz(-1.9357598) q[3];
sx q[3];
rz(-1.3273201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94806725) q[0];
sx q[0];
rz(-2.0663517) q[0];
sx q[0];
rz(1.2731592) q[0];
rz(-2.676447) q[1];
sx q[1];
rz(-1.7756614) q[1];
sx q[1];
rz(-0.21496162) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5962385) q[0];
sx q[0];
rz(-0.69850498) q[0];
sx q[0];
rz(-2.7424988) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6097492) q[2];
sx q[2];
rz(-1.1924679) q[2];
sx q[2];
rz(1.391562) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4196288) q[1];
sx q[1];
rz(-2.2441494) q[1];
sx q[1];
rz(-2.2339905) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.39542012) q[3];
sx q[3];
rz(-0.59874615) q[3];
sx q[3];
rz(1.0637763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.8421858) q[2];
sx q[2];
rz(-2.3278548) q[2];
sx q[2];
rz(2.1827533) q[2];
rz(1.1622608) q[3];
sx q[3];
rz(-1.4872888) q[3];
sx q[3];
rz(-1.6963262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6017629) q[0];
sx q[0];
rz(-1.6015263) q[0];
sx q[0];
rz(1.482561) q[0];
rz(0.70855793) q[1];
sx q[1];
rz(-2.9500912) q[1];
sx q[1];
rz(-0.74831829) q[1];
rz(-0.49025771) q[2];
sx q[2];
rz(-0.62050642) q[2];
sx q[2];
rz(-1.4668589) q[2];
rz(0.72193969) q[3];
sx q[3];
rz(-1.4293213) q[3];
sx q[3];
rz(-2.8534129) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
