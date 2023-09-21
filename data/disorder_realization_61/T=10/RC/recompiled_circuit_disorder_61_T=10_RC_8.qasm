OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7282038) q[0];
sx q[0];
rz(-2.0079186) q[0];
sx q[0];
rz(-1.5925621) q[0];
rz(-1.449861) q[1];
sx q[1];
rz(-2.4843042) q[1];
sx q[1];
rz(2.5974098) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13574164) q[0];
sx q[0];
rz(-3.0648181) q[0];
sx q[0];
rz(-0.49850054) q[0];
rz(-pi) q[1];
rz(-1.5266225) q[2];
sx q[2];
rz(-1.6720811) q[2];
sx q[2];
rz(3.0169808) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4197598) q[1];
sx q[1];
rz(-2.3082323) q[1];
sx q[1];
rz(1.014773) q[1];
rz(1.5884433) q[3];
sx q[3];
rz(-1.5994161) q[3];
sx q[3];
rz(-1.0877522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0698174) q[2];
sx q[2];
rz(-1.8775619) q[2];
sx q[2];
rz(-1.3624181) q[2];
rz(3.1133364) q[3];
sx q[3];
rz(-1.3794206) q[3];
sx q[3];
rz(0.79022592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64489275) q[0];
sx q[0];
rz(-0.70514482) q[0];
sx q[0];
rz(1.171296) q[0];
rz(0.21121875) q[1];
sx q[1];
rz(-0.44208458) q[1];
sx q[1];
rz(1.404095) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0728714) q[0];
sx q[0];
rz(-2.255548) q[0];
sx q[0];
rz(1.0147592) q[0];
x q[1];
rz(2.6589083) q[2];
sx q[2];
rz(-1.4233372) q[2];
sx q[2];
rz(2.5978616) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.055981936) q[1];
sx q[1];
rz(-1.2877052) q[1];
sx q[1];
rz(1.594747) q[1];
rz(-pi) q[2];
rz(-2.8849765) q[3];
sx q[3];
rz(-0.49391541) q[3];
sx q[3];
rz(1.6625422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.30963787) q[2];
sx q[2];
rz(-1.6654623) q[2];
sx q[2];
rz(-2.1976166) q[2];
rz(2.5850463) q[3];
sx q[3];
rz(-2.6440933) q[3];
sx q[3];
rz(-1.336162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8495162) q[0];
sx q[0];
rz(-2.3777666) q[0];
sx q[0];
rz(-1.7180432) q[0];
rz(2.3220093) q[1];
sx q[1];
rz(-0.81033605) q[1];
sx q[1];
rz(-0.56366411) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99479988) q[0];
sx q[0];
rz(-2.1437862) q[0];
sx q[0];
rz(-1.6550199) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8176739) q[2];
sx q[2];
rz(-0.62506667) q[2];
sx q[2];
rz(2.6697086) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3445661) q[1];
sx q[1];
rz(-2.1246506) q[1];
sx q[1];
rz(-0.38573854) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.54791252) q[3];
sx q[3];
rz(-0.69763819) q[3];
sx q[3];
rz(1.1988977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.50239572) q[2];
sx q[2];
rz(-1.1221308) q[2];
sx q[2];
rz(-1.0602661) q[2];
rz(3.0660196) q[3];
sx q[3];
rz(-2.2836756) q[3];
sx q[3];
rz(0.11463541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8183427) q[0];
sx q[0];
rz(-0.30775726) q[0];
sx q[0];
rz(0.19317214) q[0];
rz(1.5974143) q[1];
sx q[1];
rz(-1.9924106) q[1];
sx q[1];
rz(-2.4838122) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0021129) q[0];
sx q[0];
rz(-1.6008899) q[0];
sx q[0];
rz(-0.17480236) q[0];
rz(-pi) q[1];
rz(-2.8305956) q[2];
sx q[2];
rz(-1.9823091) q[2];
sx q[2];
rz(-2.3222773) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8448338) q[1];
sx q[1];
rz(-1.8246324) q[1];
sx q[1];
rz(2.1162507) q[1];
rz(-pi) q[2];
rz(-0.85122078) q[3];
sx q[3];
rz(-0.38098601) q[3];
sx q[3];
rz(-1.756543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0956991) q[2];
sx q[2];
rz(-0.4102439) q[2];
sx q[2];
rz(-2.0945385) q[2];
rz(-1.0632769) q[3];
sx q[3];
rz(-1.1626817) q[3];
sx q[3];
rz(1.8803546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7970153) q[0];
sx q[0];
rz(-2.8224967) q[0];
sx q[0];
rz(2.1176594) q[0];
rz(2.3539885) q[1];
sx q[1];
rz(-1.5658295) q[1];
sx q[1];
rz(-0.62686282) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.244827) q[0];
sx q[0];
rz(-1.5540431) q[0];
sx q[0];
rz(-1.5411975) q[0];
rz(-pi) q[1];
rz(-2.3847694) q[2];
sx q[2];
rz(-2.7177817) q[2];
sx q[2];
rz(-0.28048453) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.353646) q[1];
sx q[1];
rz(-1.8783356) q[1];
sx q[1];
rz(-0.31378444) q[1];
x q[2];
rz(-1.9634788) q[3];
sx q[3];
rz(-1.6851478) q[3];
sx q[3];
rz(1.8161731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.91288599) q[2];
sx q[2];
rz(-1.2330981) q[2];
sx q[2];
rz(2.1234925) q[2];
rz(2.5300238) q[3];
sx q[3];
rz(-1.3221778) q[3];
sx q[3];
rz(2.1900246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6927032) q[0];
sx q[0];
rz(-1.792181) q[0];
sx q[0];
rz(-1.1992136) q[0];
rz(1.1692283) q[1];
sx q[1];
rz(-2.0904082) q[1];
sx q[1];
rz(-2.8170524) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3015315) q[0];
sx q[0];
rz(-1.9948298) q[0];
sx q[0];
rz(-2.1246987) q[0];
x q[1];
rz(-0.80771031) q[2];
sx q[2];
rz(-1.0683904) q[2];
sx q[2];
rz(3.0903357) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.868184) q[1];
sx q[1];
rz(-1.5955828) q[1];
sx q[1];
rz(-1.0367869) q[1];
rz(-pi) q[2];
rz(0.59431521) q[3];
sx q[3];
rz(-2.6404877) q[3];
sx q[3];
rz(-1.2831812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.67363182) q[2];
sx q[2];
rz(-1.3053852) q[2];
sx q[2];
rz(1.8661873) q[2];
rz(1.0547137) q[3];
sx q[3];
rz(-2.349699) q[3];
sx q[3];
rz(-1.5462497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6773029) q[0];
sx q[0];
rz(-1.3339366) q[0];
sx q[0];
rz(1.4087079) q[0];
rz(2.6858221) q[1];
sx q[1];
rz(-0.20142889) q[1];
sx q[1];
rz(1.2021525) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31583187) q[0];
sx q[0];
rz(-1.9681265) q[0];
sx q[0];
rz(0.8855008) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1832778) q[2];
sx q[2];
rz(-0.85867184) q[2];
sx q[2];
rz(-2.9930263) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0899883) q[1];
sx q[1];
rz(-3.0172536) q[1];
sx q[1];
rz(-0.73595562) q[1];
rz(-pi) q[2];
rz(2.8771411) q[3];
sx q[3];
rz(-1.7489986) q[3];
sx q[3];
rz(1.4101392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.000164) q[2];
sx q[2];
rz(-1.8464073) q[2];
sx q[2];
rz(-0.84623519) q[2];
rz(-2.6464461) q[3];
sx q[3];
rz(-0.64619243) q[3];
sx q[3];
rz(-1.600986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0705567) q[0];
sx q[0];
rz(-1.2905916) q[0];
sx q[0];
rz(1.113168) q[0];
rz(2.44599) q[1];
sx q[1];
rz(-0.38989392) q[1];
sx q[1];
rz(-1.5323458) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1217864) q[0];
sx q[0];
rz(-1.1018503) q[0];
sx q[0];
rz(1.8437587) q[0];
x q[1];
rz(1.1148909) q[2];
sx q[2];
rz(-2.0156983) q[2];
sx q[2];
rz(2.0632495) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.53837) q[1];
sx q[1];
rz(-2.6909628) q[1];
sx q[1];
rz(1.6166812) q[1];
rz(-pi) q[2];
rz(-1.6892151) q[3];
sx q[3];
rz(-2.1695671) q[3];
sx q[3];
rz(1.0373725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9939076) q[2];
sx q[2];
rz(-0.2299749) q[2];
sx q[2];
rz(0.85419401) q[2];
rz(1.9130075) q[3];
sx q[3];
rz(-1.5766141) q[3];
sx q[3];
rz(1.6555697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5937186) q[0];
sx q[0];
rz(-0.49867189) q[0];
sx q[0];
rz(-0.6643995) q[0];
rz(-1.1876855) q[1];
sx q[1];
rz(-0.70078754) q[1];
sx q[1];
rz(1.857035) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5212095) q[0];
sx q[0];
rz(-1.2854344) q[0];
sx q[0];
rz(-0.83831212) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.2662439) q[2];
sx q[2];
rz(-1.2028482) q[2];
sx q[2];
rz(3.0014696) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7582693) q[1];
sx q[1];
rz(-2.3136534) q[1];
sx q[1];
rz(0.62288706) q[1];
x q[2];
rz(-1.5434389) q[3];
sx q[3];
rz(-1.9606855) q[3];
sx q[3];
rz(-2.011812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5294042) q[2];
sx q[2];
rz(-2.2470784) q[2];
sx q[2];
rz(-2.5908296) q[2];
rz(-2.4380056) q[3];
sx q[3];
rz(-2.1313322) q[3];
sx q[3];
rz(1.5065058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4187014) q[0];
sx q[0];
rz(-0.35690618) q[0];
sx q[0];
rz(-0.044145949) q[0];
rz(-1.528953) q[1];
sx q[1];
rz(-1.9020558) q[1];
sx q[1];
rz(0.77967656) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6990307) q[0];
sx q[0];
rz(-0.5578707) q[0];
sx q[0];
rz(2.872422) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.81515628) q[2];
sx q[2];
rz(-1.2218468) q[2];
sx q[2];
rz(2.7100035) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.979094) q[1];
sx q[1];
rz(-1.7531803) q[1];
sx q[1];
rz(2.6136293) q[1];
rz(-pi) q[2];
rz(2.3603504) q[3];
sx q[3];
rz(-1.9133854) q[3];
sx q[3];
rz(-0.70171802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6471275) q[2];
sx q[2];
rz(-2.4202042) q[2];
sx q[2];
rz(1.5987827) q[2];
rz(-2.4370082) q[3];
sx q[3];
rz(-1.6274118) q[3];
sx q[3];
rz(-0.012044756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4298532) q[0];
sx q[0];
rz(-1.585351) q[0];
sx q[0];
rz(2.4899695) q[0];
rz(1.7977057) q[1];
sx q[1];
rz(-1.4812891) q[1];
sx q[1];
rz(-0.67566009) q[1];
rz(-1.1336397) q[2];
sx q[2];
rz(-1.8532248) q[2];
sx q[2];
rz(2.1291477) q[2];
rz(2.5034954) q[3];
sx q[3];
rz(-1.6774072) q[3];
sx q[3];
rz(0.93117136) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];