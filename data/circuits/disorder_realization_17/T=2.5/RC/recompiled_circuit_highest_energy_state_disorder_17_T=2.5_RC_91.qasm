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
rz(1.0139581) q[0];
sx q[0];
rz(-1.6887709) q[0];
sx q[0];
rz(-0.59803522) q[0];
rz(1.1818089) q[1];
sx q[1];
rz(-0.66507566) q[1];
sx q[1];
rz(-2.9887078) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0273697) q[0];
sx q[0];
rz(-1.4589757) q[0];
sx q[0];
rz(-3.0614475) q[0];
x q[1];
rz(3.0316584) q[2];
sx q[2];
rz(-2.6510932) q[2];
sx q[2];
rz(-2.5562499) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.26273706) q[1];
sx q[1];
rz(-1.3595743) q[1];
sx q[1];
rz(-3.0131574) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0956354) q[3];
sx q[3];
rz(-1.5737688) q[3];
sx q[3];
rz(-1.3770905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.87024706) q[2];
sx q[2];
rz(-2.953244) q[2];
sx q[2];
rz(1.9625473) q[2];
rz(-0.41006586) q[3];
sx q[3];
rz(-2.0867917) q[3];
sx q[3];
rz(-2.2470391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9533185) q[0];
sx q[0];
rz(-0.66903791) q[0];
sx q[0];
rz(2.8642995) q[0];
rz(0.20031985) q[1];
sx q[1];
rz(-2.2453313) q[1];
sx q[1];
rz(-0.083981363) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6407627) q[0];
sx q[0];
rz(-2.6272054) q[0];
sx q[0];
rz(-1.1577093) q[0];
rz(-pi) q[1];
rz(2.9756184) q[2];
sx q[2];
rz(-0.73455185) q[2];
sx q[2];
rz(-2.4935472) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9878432) q[1];
sx q[1];
rz(-2.0122408) q[1];
sx q[1];
rz(0.50905871) q[1];
rz(-1.3091607) q[3];
sx q[3];
rz(-2.0896834) q[3];
sx q[3];
rz(2.6760657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1818485) q[2];
sx q[2];
rz(-2.4740348) q[2];
sx q[2];
rz(-2.3954771) q[2];
rz(-0.6212081) q[3];
sx q[3];
rz(-0.64380232) q[3];
sx q[3];
rz(1.752689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-2.8212432) q[0];
sx q[0];
rz(-0.72750434) q[0];
sx q[0];
rz(1.2028836) q[0];
rz(-2.7299643) q[1];
sx q[1];
rz(-1.2806712) q[1];
sx q[1];
rz(-2.0278377) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4899556) q[0];
sx q[0];
rz(-2.2764766) q[0];
sx q[0];
rz(2.5791427) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.95543785) q[2];
sx q[2];
rz(-1.9410543) q[2];
sx q[2];
rz(0.40981612) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3948156) q[1];
sx q[1];
rz(-1.4377497) q[1];
sx q[1];
rz(-2.2486671) q[1];
x q[2];
rz(2.2334072) q[3];
sx q[3];
rz(-0.81530967) q[3];
sx q[3];
rz(2.1774928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0456298) q[2];
sx q[2];
rz(-0.7926422) q[2];
sx q[2];
rz(-1.1620129) q[2];
rz(-2.3368733) q[3];
sx q[3];
rz(-1.8685124) q[3];
sx q[3];
rz(1.8603676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3941037) q[0];
sx q[0];
rz(-1.2825613) q[0];
sx q[0];
rz(1.5149186) q[0];
rz(2.6331242) q[1];
sx q[1];
rz(-0.6503121) q[1];
sx q[1];
rz(-1.633684) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28428299) q[0];
sx q[0];
rz(-1.4197667) q[0];
sx q[0];
rz(3.1394464) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1437683) q[2];
sx q[2];
rz(-0.25616872) q[2];
sx q[2];
rz(-1.0333956) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0421988) q[1];
sx q[1];
rz(-1.458199) q[1];
sx q[1];
rz(-1.8774752) q[1];
x q[2];
rz(1.0359457) q[3];
sx q[3];
rz(-0.7300762) q[3];
sx q[3];
rz(1.4393782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7278829) q[2];
sx q[2];
rz(-1.1943613) q[2];
sx q[2];
rz(0.74753648) q[2];
rz(-1.9109292) q[3];
sx q[3];
rz(-2.3321798) q[3];
sx q[3];
rz(-2.6443853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.704945) q[0];
sx q[0];
rz(-2.0132988) q[0];
sx q[0];
rz(-1.6823912) q[0];
rz(2.7875426) q[1];
sx q[1];
rz(-0.67146462) q[1];
sx q[1];
rz(-0.57463247) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.282772) q[0];
sx q[0];
rz(-3.0631815) q[0];
sx q[0];
rz(0.17228244) q[0];
x q[1];
rz(-1.2990405) q[2];
sx q[2];
rz(-0.82052204) q[2];
sx q[2];
rz(-2.3472903) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1916442) q[1];
sx q[1];
rz(-0.23298082) q[1];
sx q[1];
rz(2.9162372) q[1];
rz(-pi) q[2];
rz(2.0784723) q[3];
sx q[3];
rz(-2.9083544) q[3];
sx q[3];
rz(2.4539029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6358801) q[2];
sx q[2];
rz(-1.1735703) q[2];
sx q[2];
rz(-2.9393348) q[2];
rz(-2.108719) q[3];
sx q[3];
rz(-2.5346916) q[3];
sx q[3];
rz(-3.0752944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1101087) q[0];
sx q[0];
rz(-0.38723543) q[0];
sx q[0];
rz(0.37145823) q[0];
rz(-2.8455041) q[1];
sx q[1];
rz(-1.5541872) q[1];
sx q[1];
rz(2.9782226) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4811081) q[0];
sx q[0];
rz(-2.6077932) q[0];
sx q[0];
rz(0.36719881) q[0];
rz(-pi) q[1];
rz(0.67739112) q[2];
sx q[2];
rz(-2.4755726) q[2];
sx q[2];
rz(2.1474554) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5886053) q[1];
sx q[1];
rz(-1.7514015) q[1];
sx q[1];
rz(-0.44579472) q[1];
rz(-pi) q[2];
rz(-1.7512239) q[3];
sx q[3];
rz(-2.437915) q[3];
sx q[3];
rz(-2.0247519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0614193) q[2];
sx q[2];
rz(-0.27013186) q[2];
sx q[2];
rz(0.6529676) q[2];
rz(-0.4717007) q[3];
sx q[3];
rz(-0.74840778) q[3];
sx q[3];
rz(-0.72699839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2111874) q[0];
sx q[0];
rz(-1.0791595) q[0];
sx q[0];
rz(2.1042673) q[0];
rz(0.51517454) q[1];
sx q[1];
rz(-1.2778792) q[1];
sx q[1];
rz(-1.8151262) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7221092) q[0];
sx q[0];
rz(-0.67841716) q[0];
sx q[0];
rz(-2.2161311) q[0];
rz(-pi) q[1];
x q[1];
rz(0.35197978) q[2];
sx q[2];
rz(-1.5328141) q[2];
sx q[2];
rz(0.69881907) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.675589) q[1];
sx q[1];
rz(-1.4741305) q[1];
sx q[1];
rz(-2.3748891) q[1];
x q[2];
rz(-1.8765728) q[3];
sx q[3];
rz(-1.5103087) q[3];
sx q[3];
rz(-2.4779589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4482164) q[2];
sx q[2];
rz(-2.6770112) q[2];
sx q[2];
rz(2.8153937) q[2];
rz(2.163573) q[3];
sx q[3];
rz(-1.2468485) q[3];
sx q[3];
rz(3.0514362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6982894) q[0];
sx q[0];
rz(-0.25743085) q[0];
sx q[0];
rz(0.28939104) q[0];
rz(3.1293213) q[1];
sx q[1];
rz(-2.8616276) q[1];
sx q[1];
rz(-1.6562921) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12690954) q[0];
sx q[0];
rz(-2.2193847) q[0];
sx q[0];
rz(-2.4832583) q[0];
rz(-pi) q[1];
rz(1.6248762) q[2];
sx q[2];
rz(-1.7940494) q[2];
sx q[2];
rz(0.28424451) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.43013182) q[1];
sx q[1];
rz(-1.3688812) q[1];
sx q[1];
rz(1.422775) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1505707) q[3];
sx q[3];
rz(-0.54383792) q[3];
sx q[3];
rz(-0.71292927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5926008) q[2];
sx q[2];
rz(-1.5697378) q[2];
sx q[2];
rz(-0.17295095) q[2];
rz(-1.7671827) q[3];
sx q[3];
rz(-1.3783312) q[3];
sx q[3];
rz(-1.6636728) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7462815) q[0];
sx q[0];
rz(-2.6938541) q[0];
sx q[0];
rz(2.3315499) q[0];
rz(0.41656247) q[1];
sx q[1];
rz(-2.3655393) q[1];
sx q[1];
rz(2.7966444) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5589011) q[0];
sx q[0];
rz(-1.693629) q[0];
sx q[0];
rz(2.4305811) q[0];
x q[1];
rz(-1.2917963) q[2];
sx q[2];
rz(-1.7300182) q[2];
sx q[2];
rz(2.4375242) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4629412) q[1];
sx q[1];
rz(-0.57616568) q[1];
sx q[1];
rz(2.1919996) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3446001) q[3];
sx q[3];
rz(-0.16599645) q[3];
sx q[3];
rz(0.57193631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.42744669) q[2];
sx q[2];
rz(-2.2553208) q[2];
sx q[2];
rz(0.084058849) q[2];
rz(-0.90703026) q[3];
sx q[3];
rz(-1.7232938) q[3];
sx q[3];
rz(-2.1840054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8091938) q[0];
sx q[0];
rz(-3.0539303) q[0];
sx q[0];
rz(-2.8293389) q[0];
rz(-2.9088083) q[1];
sx q[1];
rz(-0.39118958) q[1];
sx q[1];
rz(1.8544474) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2256266) q[0];
sx q[0];
rz(-0.73333987) q[0];
sx q[0];
rz(0.49247964) q[0];
x q[1];
rz(1.5239117) q[2];
sx q[2];
rz(-2.2094036) q[2];
sx q[2];
rz(0.70568161) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6010661) q[1];
sx q[1];
rz(-2.6208715) q[1];
sx q[1];
rz(-2.0784318) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.62508836) q[3];
sx q[3];
rz(-1.8755138) q[3];
sx q[3];
rz(-2.8947322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9707668) q[2];
sx q[2];
rz(-1.2831186) q[2];
sx q[2];
rz(1.4975366) q[2];
rz(1.0981285) q[3];
sx q[3];
rz(-0.69881717) q[3];
sx q[3];
rz(0.48925492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.004414) q[0];
sx q[0];
rz(-1.7902086) q[0];
sx q[0];
rz(-0.93850346) q[0];
rz(1.3068403) q[1];
sx q[1];
rz(-2.2521781) q[1];
sx q[1];
rz(2.3440012) q[1];
rz(-2.9970279) q[2];
sx q[2];
rz(-2.2546386) q[2];
sx q[2];
rz(-2.4972514) q[2];
rz(-0.90436892) q[3];
sx q[3];
rz(-1.6828129) q[3];
sx q[3];
rz(-3.070698) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
