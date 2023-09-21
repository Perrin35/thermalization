OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9053469) q[0];
sx q[0];
rz(5.5570931) q[0];
sx q[0];
rz(9.2232016) q[0];
rz(-2.6456614) q[1];
sx q[1];
rz(-2.6013241) q[1];
sx q[1];
rz(0.93710605) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9286081) q[0];
sx q[0];
rz(-2.1323418) q[0];
sx q[0];
rz(-2.0748078) q[0];
x q[1];
rz(-2.5295528) q[2];
sx q[2];
rz(-1.1571552) q[2];
sx q[2];
rz(0.94758247) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5697437) q[1];
sx q[1];
rz(-1.3880001) q[1];
sx q[1];
rz(0.079449541) q[1];
rz(-1.8398383) q[3];
sx q[3];
rz(-1.4770513) q[3];
sx q[3];
rz(-2.4274488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9822838) q[2];
sx q[2];
rz(-0.26580492) q[2];
sx q[2];
rz(0.086159555) q[2];
rz(0.75749767) q[3];
sx q[3];
rz(-0.75452724) q[3];
sx q[3];
rz(1.0936201) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5646097) q[0];
sx q[0];
rz(-1.4819205) q[0];
sx q[0];
rz(1.8151059) q[0];
rz(-1.2558698) q[1];
sx q[1];
rz(-1.5763667) q[1];
sx q[1];
rz(0.27145162) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46389929) q[0];
sx q[0];
rz(-1.1743744) q[0];
sx q[0];
rz(-2.1167614) q[0];
rz(-pi) q[1];
rz(1.2618622) q[2];
sx q[2];
rz(-2.0663107) q[2];
sx q[2];
rz(2.3625284) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2084864) q[1];
sx q[1];
rz(-2.3098364) q[1];
sx q[1];
rz(-1.0528802) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3321757) q[3];
sx q[3];
rz(-1.8339694) q[3];
sx q[3];
rz(-1.1337048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0945956) q[2];
sx q[2];
rz(-1.3826028) q[2];
sx q[2];
rz(0.57717741) q[2];
rz(0.92352891) q[3];
sx q[3];
rz(-0.90463343) q[3];
sx q[3];
rz(1.4097759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.24401027) q[0];
sx q[0];
rz(-1.0616466) q[0];
sx q[0];
rz(0.14257167) q[0];
rz(-1.7890731) q[1];
sx q[1];
rz(-1.0357772) q[1];
sx q[1];
rz(-0.20908633) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36172983) q[0];
sx q[0];
rz(-1.0233101) q[0];
sx q[0];
rz(-2.2453383) q[0];
rz(-pi) q[1];
rz(-1.9833343) q[2];
sx q[2];
rz(-0.87450714) q[2];
sx q[2];
rz(-1.0227026) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0774539) q[1];
sx q[1];
rz(-1.9874548) q[1];
sx q[1];
rz(2.0616848) q[1];
x q[2];
rz(1.7861869) q[3];
sx q[3];
rz(-1.2743909) q[3];
sx q[3];
rz(-2.4733558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3645939) q[2];
sx q[2];
rz(-0.89407388) q[2];
sx q[2];
rz(0.40412942) q[2];
rz(-1.8557619) q[3];
sx q[3];
rz(-2.0127681) q[3];
sx q[3];
rz(-0.16734853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4362815) q[0];
sx q[0];
rz(-0.10665882) q[0];
sx q[0];
rz(-2.1787815) q[0];
rz(-0.46936938) q[1];
sx q[1];
rz(-0.58987394) q[1];
sx q[1];
rz(0.00096360047) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22165933) q[0];
sx q[0];
rz(-0.98826212) q[0];
sx q[0];
rz(-0.68908738) q[0];
rz(3.0175866) q[2];
sx q[2];
rz(-1.5702855) q[2];
sx q[2];
rz(-0.83204568) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.094593781) q[1];
sx q[1];
rz(-1.440289) q[1];
sx q[1];
rz(0.99884896) q[1];
rz(-pi) q[2];
rz(-1.4812874) q[3];
sx q[3];
rz(-1.0032017) q[3];
sx q[3];
rz(-1.2168509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0756388) q[2];
sx q[2];
rz(-1.5116189) q[2];
sx q[2];
rz(-0.11165079) q[2];
rz(-2.3305437) q[3];
sx q[3];
rz(-2.6871197) q[3];
sx q[3];
rz(3.1276935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1902996) q[0];
sx q[0];
rz(-0.90068978) q[0];
sx q[0];
rz(2.7850889) q[0];
rz(0.50645343) q[1];
sx q[1];
rz(-1.3566147) q[1];
sx q[1];
rz(0.26062632) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60359943) q[0];
sx q[0];
rz(-1.7126181) q[0];
sx q[0];
rz(1.5285368) q[0];
rz(-pi) q[1];
rz(1.8243276) q[2];
sx q[2];
rz(-1.8142482) q[2];
sx q[2];
rz(-1.4075116) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4958447) q[1];
sx q[1];
rz(-2.145088) q[1];
sx q[1];
rz(1.5498284) q[1];
rz(1.4277677) q[3];
sx q[3];
rz(-1.7149394) q[3];
sx q[3];
rz(0.10728697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9115209) q[2];
sx q[2];
rz(-1.4804966) q[2];
sx q[2];
rz(-2.9122706) q[2];
rz(2.5991332) q[3];
sx q[3];
rz(-2.8312603) q[3];
sx q[3];
rz(2.1054161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3689573) q[0];
sx q[0];
rz(-2.2326523) q[0];
sx q[0];
rz(-1.4916346) q[0];
rz(-1.0391327) q[1];
sx q[1];
rz(-1.2960478) q[1];
sx q[1];
rz(-1.7274436) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3154253) q[0];
sx q[0];
rz(-2.8624479) q[0];
sx q[0];
rz(-2.9690353) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.44712375) q[2];
sx q[2];
rz(-0.56338718) q[2];
sx q[2];
rz(2.2396357) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.27281877) q[1];
sx q[1];
rz(-2.0588074) q[1];
sx q[1];
rz(-0.52656071) q[1];
rz(-pi) q[2];
x q[2];
rz(3.066091) q[3];
sx q[3];
rz(-1.7864979) q[3];
sx q[3];
rz(0.84738934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1386537) q[2];
sx q[2];
rz(-2.5230375) q[2];
sx q[2];
rz(3.1138528) q[2];
rz(2.6489143) q[3];
sx q[3];
rz(-1.6511107) q[3];
sx q[3];
rz(1.3940575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3843) q[0];
sx q[0];
rz(-1.5526271) q[0];
sx q[0];
rz(-0.12167715) q[0];
rz(-1.9901468) q[1];
sx q[1];
rz(-2.6897488) q[1];
sx q[1];
rz(-2.7391403) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2936195) q[0];
sx q[0];
rz(-1.3012039) q[0];
sx q[0];
rz(-0.86975354) q[0];
x q[1];
rz(1.5624814) q[2];
sx q[2];
rz(-1.6168211) q[2];
sx q[2];
rz(-0.41551057) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.23558815) q[1];
sx q[1];
rz(-1.2885805) q[1];
sx q[1];
rz(1.7547592) q[1];
rz(-pi) q[2];
rz(-1.6472858) q[3];
sx q[3];
rz(-0.34919958) q[3];
sx q[3];
rz(2.4793712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1272614) q[2];
sx q[2];
rz(-1.971259) q[2];
sx q[2];
rz(0.86223117) q[2];
rz(2.6640653) q[3];
sx q[3];
rz(-1.9600441) q[3];
sx q[3];
rz(-0.91807085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7858793) q[0];
sx q[0];
rz(-0.15545758) q[0];
sx q[0];
rz(-2.4086337) q[0];
rz(0.14006242) q[1];
sx q[1];
rz(-2.1439794) q[1];
sx q[1];
rz(2.1070811) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6013992) q[0];
sx q[0];
rz(-1.7978151) q[0];
sx q[0];
rz(2.9674203) q[0];
x q[1];
rz(-1.6413692) q[2];
sx q[2];
rz(-1.4683873) q[2];
sx q[2];
rz(0.18825738) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.95755277) q[1];
sx q[1];
rz(-1.5337481) q[1];
sx q[1];
rz(2.1368105) q[1];
rz(-pi) q[2];
rz(0.23969527) q[3];
sx q[3];
rz(-0.42907676) q[3];
sx q[3];
rz(-2.808188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.90116477) q[2];
sx q[2];
rz(-1.1952091) q[2];
sx q[2];
rz(-0.77587664) q[2];
rz(2.4173229) q[3];
sx q[3];
rz(-2.3359719) q[3];
sx q[3];
rz(1.7075214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4416606) q[0];
sx q[0];
rz(-0.76403809) q[0];
sx q[0];
rz(0.016816703) q[0];
rz(3.1230714) q[1];
sx q[1];
rz(-1.3341981) q[1];
sx q[1];
rz(2.3628078) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.461207) q[0];
sx q[0];
rz(-1.2347504) q[0];
sx q[0];
rz(1.3075605) q[0];
rz(-pi) q[1];
rz(-1.0672827) q[2];
sx q[2];
rz(-0.52769606) q[2];
sx q[2];
rz(2.6380981) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3719337) q[1];
sx q[1];
rz(-1.2049335) q[1];
sx q[1];
rz(-2.8918173) q[1];
x q[2];
rz(1.7917463) q[3];
sx q[3];
rz(-1.5077935) q[3];
sx q[3];
rz(2.8642879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4497711) q[2];
sx q[2];
rz(-1.8372953) q[2];
sx q[2];
rz(2.4251535) q[2];
rz(-1.6843494) q[3];
sx q[3];
rz(-1.4102035) q[3];
sx q[3];
rz(1.8523857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0338106) q[0];
sx q[0];
rz(-2.6066715) q[0];
sx q[0];
rz(2.9901436) q[0];
rz(-2.9653213) q[1];
sx q[1];
rz(-1.9468032) q[1];
sx q[1];
rz(-0.72296468) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13847362) q[0];
sx q[0];
rz(-1.4644633) q[0];
sx q[0];
rz(1.4195725) q[0];
rz(1.1935913) q[2];
sx q[2];
rz(-1.4133487) q[2];
sx q[2];
rz(2.8709656) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.55262676) q[1];
sx q[1];
rz(-0.50682658) q[1];
sx q[1];
rz(1.6848906) q[1];
x q[2];
rz(-0.39977269) q[3];
sx q[3];
rz(-0.82953605) q[3];
sx q[3];
rz(-1.6403891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.67939776) q[2];
sx q[2];
rz(-1.858819) q[2];
sx q[2];
rz(2.6160713) q[2];
rz(2.8578791) q[3];
sx q[3];
rz(-1.107639) q[3];
sx q[3];
rz(2.7450558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5491966) q[0];
sx q[0];
rz(-1.1239197) q[0];
sx q[0];
rz(-0.32604937) q[0];
rz(1.0992959) q[1];
sx q[1];
rz(-0.14930832) q[1];
sx q[1];
rz(2.2717448) q[1];
rz(-2.6635086) q[2];
sx q[2];
rz(-2.1838084) q[2];
sx q[2];
rz(2.9391391) q[2];
rz(1.7952193) q[3];
sx q[3];
rz(-1.2524458) q[3];
sx q[3];
rz(2.8129775) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];