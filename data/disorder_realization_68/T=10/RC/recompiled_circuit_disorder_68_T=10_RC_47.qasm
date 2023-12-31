OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.39419898) q[0];
sx q[0];
rz(-0.49180254) q[0];
sx q[0];
rz(-2.9536182) q[0];
rz(-1.1176874) q[1];
sx q[1];
rz(-1.517065) q[1];
sx q[1];
rz(2.7741073) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1023094) q[0];
sx q[0];
rz(-1.4674205) q[0];
sx q[0];
rz(-2.1084059) q[0];
rz(-pi) q[1];
rz(0.1349749) q[2];
sx q[2];
rz(-1.0834603) q[2];
sx q[2];
rz(-0.48263532) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6403113) q[1];
sx q[1];
rz(-2.8456563) q[1];
sx q[1];
rz(-0.96247767) q[1];
rz(-pi) q[2];
rz(0.46842694) q[3];
sx q[3];
rz(-2.763063) q[3];
sx q[3];
rz(2.1551876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.964103) q[2];
sx q[2];
rz(-0.51012817) q[2];
sx q[2];
rz(-0.5509848) q[2];
rz(-1.8356813) q[3];
sx q[3];
rz(-1.6492313) q[3];
sx q[3];
rz(1.8252385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
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
rz(0.47857639) q[0];
sx q[0];
rz(-1.0000279) q[0];
sx q[0];
rz(-0.4719032) q[0];
rz(2.7117803) q[1];
sx q[1];
rz(-1.2496354) q[1];
sx q[1];
rz(0.93634161) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3204526) q[0];
sx q[0];
rz(-1.8006514) q[0];
sx q[0];
rz(-1.6605404) q[0];
rz(-pi) q[1];
rz(1.9905375) q[2];
sx q[2];
rz(-2.310576) q[2];
sx q[2];
rz(-1.3295528) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4221103) q[1];
sx q[1];
rz(-1.6750095) q[1];
sx q[1];
rz(-2.4476493) q[1];
rz(-0.02336054) q[3];
sx q[3];
rz(-2.0733895) q[3];
sx q[3];
rz(2.392829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.77461809) q[2];
sx q[2];
rz(-0.32704157) q[2];
sx q[2];
rz(2.7152087) q[2];
rz(1.2373699) q[3];
sx q[3];
rz(-2.5137413) q[3];
sx q[3];
rz(0.0330851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8957829) q[0];
sx q[0];
rz(-1.824546) q[0];
sx q[0];
rz(0.93908969) q[0];
rz(-0.89871961) q[1];
sx q[1];
rz(-0.4788613) q[1];
sx q[1];
rz(2.5476707) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31375162) q[0];
sx q[0];
rz(-1.272164) q[0];
sx q[0];
rz(-1.2350425) q[0];
rz(-pi) q[1];
rz(2.2592696) q[2];
sx q[2];
rz(-1.6154628) q[2];
sx q[2];
rz(0.67827144) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8201901) q[1];
sx q[1];
rz(-0.47649511) q[1];
sx q[1];
rz(2.0501775) q[1];
x q[2];
rz(1.9150312) q[3];
sx q[3];
rz(-2.0153181) q[3];
sx q[3];
rz(-1.3821186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.64017355) q[2];
sx q[2];
rz(-2.4121425) q[2];
sx q[2];
rz(1.4397941) q[2];
rz(2.7539608) q[3];
sx q[3];
rz(-1.6250316) q[3];
sx q[3];
rz(0.38813996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7664292) q[0];
sx q[0];
rz(-1.5513993) q[0];
sx q[0];
rz(2.638812) q[0];
rz(0.76820961) q[1];
sx q[1];
rz(-2.6380824) q[1];
sx q[1];
rz(2.3847413) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6605646) q[0];
sx q[0];
rz(-1.2111944) q[0];
sx q[0];
rz(-1.9268376) q[0];
x q[1];
rz(-1.5966162) q[2];
sx q[2];
rz(-2.6757247) q[2];
sx q[2];
rz(2.4398838) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.56843578) q[1];
sx q[1];
rz(-0.70306289) q[1];
sx q[1];
rz(-0.99848024) q[1];
x q[2];
rz(-2.5960856) q[3];
sx q[3];
rz(-2.3295998) q[3];
sx q[3];
rz(-0.10520392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.42671529) q[2];
sx q[2];
rz(-1.9116414) q[2];
sx q[2];
rz(1.654401) q[2];
rz(0.58250827) q[3];
sx q[3];
rz(-2.0472066) q[3];
sx q[3];
rz(-0.55707651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8476167) q[0];
sx q[0];
rz(-1.0852381) q[0];
sx q[0];
rz(2.3838682) q[0];
rz(1.2879397) q[1];
sx q[1];
rz(-2.2133591) q[1];
sx q[1];
rz(1.0505189) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8331063) q[0];
sx q[0];
rz(-1.8252488) q[0];
sx q[0];
rz(1.2480877) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1804785) q[2];
sx q[2];
rz(-0.69176199) q[2];
sx q[2];
rz(1.320653) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1615636) q[1];
sx q[1];
rz(-1.0161576) q[1];
sx q[1];
rz(2.5142923) q[1];
rz(-2.3349808) q[3];
sx q[3];
rz(-1.9730554) q[3];
sx q[3];
rz(0.81392399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.918255) q[2];
sx q[2];
rz(-2.788322) q[2];
sx q[2];
rz(-0.63344947) q[2];
rz(-1.194362) q[3];
sx q[3];
rz(-1.6703689) q[3];
sx q[3];
rz(0.71715322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69960064) q[0];
sx q[0];
rz(-0.49015912) q[0];
sx q[0];
rz(-0.25318405) q[0];
rz(1.5340012) q[1];
sx q[1];
rz(-1.7065159) q[1];
sx q[1];
rz(1.4621428) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7397241) q[0];
sx q[0];
rz(-1.3603633) q[0];
sx q[0];
rz(1.6519288) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2199549) q[2];
sx q[2];
rz(-1.3772794) q[2];
sx q[2];
rz(2.9606539) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.05789214) q[1];
sx q[1];
rz(-1.9533227) q[1];
sx q[1];
rz(-1.7544569) q[1];
x q[2];
rz(-1.9147647) q[3];
sx q[3];
rz(-1.2729537) q[3];
sx q[3];
rz(-0.96836585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9399461) q[2];
sx q[2];
rz(-2.3942409) q[2];
sx q[2];
rz(-0.80491006) q[2];
rz(1.1770052) q[3];
sx q[3];
rz(-1.0369119) q[3];
sx q[3];
rz(2.0578407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8686304) q[0];
sx q[0];
rz(-1.0694163) q[0];
sx q[0];
rz(-0.92765635) q[0];
rz(2.1169128) q[1];
sx q[1];
rz(-1.6352446) q[1];
sx q[1];
rz(1.0120846) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1597848) q[0];
sx q[0];
rz(-0.73043981) q[0];
sx q[0];
rz(-0.83321379) q[0];
x q[1];
rz(0.32832844) q[2];
sx q[2];
rz(-0.40847455) q[2];
sx q[2];
rz(0.35818737) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.183179) q[1];
sx q[1];
rz(-1.603754) q[1];
sx q[1];
rz(-2.5977913) q[1];
rz(0.55862553) q[3];
sx q[3];
rz(-1.6405676) q[3];
sx q[3];
rz(0.39486265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6107789) q[2];
sx q[2];
rz(-1.4799708) q[2];
sx q[2];
rz(0.42993316) q[2];
rz(2.1271465) q[3];
sx q[3];
rz(-0.40922624) q[3];
sx q[3];
rz(0.53340069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6124509) q[0];
sx q[0];
rz(-2.2021459) q[0];
sx q[0];
rz(-2.9274143) q[0];
rz(-1.051349) q[1];
sx q[1];
rz(-2.9290757) q[1];
sx q[1];
rz(-2.8578551) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8026233) q[0];
sx q[0];
rz(-1.5366652) q[0];
sx q[0];
rz(2.8404833) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2405422) q[2];
sx q[2];
rz(-1.3888161) q[2];
sx q[2];
rz(2.6957126) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2092065) q[1];
sx q[1];
rz(-1.7665518) q[1];
sx q[1];
rz(2.6372361) q[1];
x q[2];
rz(2.22105) q[3];
sx q[3];
rz(-1.0245819) q[3];
sx q[3];
rz(1.9073245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.45067898) q[2];
sx q[2];
rz(-2.6288855) q[2];
sx q[2];
rz(1.696375) q[2];
rz(1.5444267) q[3];
sx q[3];
rz(-1.4404567) q[3];
sx q[3];
rz(-2.8022695) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63672367) q[0];
sx q[0];
rz(-3.0915785) q[0];
sx q[0];
rz(3.0723363) q[0];
rz(-1.4878558) q[1];
sx q[1];
rz(-1.8585049) q[1];
sx q[1];
rz(1.5725296) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42800316) q[0];
sx q[0];
rz(-2.028095) q[0];
sx q[0];
rz(-0.86488117) q[0];
x q[1];
rz(-2.9022129) q[2];
sx q[2];
rz(-0.33472543) q[2];
sx q[2];
rz(2.7811546) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6988277) q[1];
sx q[1];
rz(-0.24833939) q[1];
sx q[1];
rz(2.793924) q[1];
rz(1.8435555) q[3];
sx q[3];
rz(-2.4747304) q[3];
sx q[3];
rz(0.71344261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1853603) q[2];
sx q[2];
rz(-2.9113443) q[2];
sx q[2];
rz(3.0017079) q[2];
rz(-2.774003) q[3];
sx q[3];
rz(-1.1871754) q[3];
sx q[3];
rz(-2.1504413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96520987) q[0];
sx q[0];
rz(-2.7503224) q[0];
sx q[0];
rz(-2.4998253) q[0];
rz(1.2311252) q[1];
sx q[1];
rz(-1.9893913) q[1];
sx q[1];
rz(0.26783255) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5689241) q[0];
sx q[0];
rz(-2.2974282) q[0];
sx q[0];
rz(0.6092682) q[0];
rz(-0.12561663) q[2];
sx q[2];
rz(-1.7526502) q[2];
sx q[2];
rz(-0.4609209) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6707582) q[1];
sx q[1];
rz(-1.8948312) q[1];
sx q[1];
rz(1.2653989) q[1];
rz(-2.0249428) q[3];
sx q[3];
rz(-0.59026679) q[3];
sx q[3];
rz(0.86436194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3315167) q[2];
sx q[2];
rz(-1.6167567) q[2];
sx q[2];
rz(1.6646741) q[2];
rz(2.8752575) q[3];
sx q[3];
rz(-0.24644066) q[3];
sx q[3];
rz(2.5951071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1289566) q[0];
sx q[0];
rz(-2.2205882) q[0];
sx q[0];
rz(0.90482774) q[0];
rz(-2.3616882) q[1];
sx q[1];
rz(-2.654568) q[1];
sx q[1];
rz(1.754896) q[1];
rz(0.71074625) q[2];
sx q[2];
rz(-1.1087316) q[2];
sx q[2];
rz(1.0089594) q[2];
rz(-0.6775425) q[3];
sx q[3];
rz(-1.7384221) q[3];
sx q[3];
rz(1.3330028) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
