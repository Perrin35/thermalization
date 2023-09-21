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
rz(2.6497901) q[0];
sx q[0];
rz(9.2368035) q[0];
rz(2.0239053) q[1];
sx q[1];
rz(-1.6245276) q[1];
sx q[1];
rz(-2.7741073) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6715235) q[0];
sx q[0];
rz(-1.0363665) q[0];
sx q[0];
rz(0.1202017) q[0];
rz(-1.0796702) q[2];
sx q[2];
rz(-1.4516146) q[2];
sx q[2];
rz(2.116938) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.130598) q[1];
sx q[1];
rz(-1.8124609) q[1];
sx q[1];
rz(-0.17250891) q[1];
rz(-1.7484619) q[3];
sx q[3];
rz(-1.2347617) q[3];
sx q[3];
rz(0.48776585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1774896) q[2];
sx q[2];
rz(-2.6314645) q[2];
sx q[2];
rz(0.5509848) q[2];
rz(-1.3059113) q[3];
sx q[3];
rz(-1.6492313) q[3];
sx q[3];
rz(-1.8252385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47857639) q[0];
sx q[0];
rz(-2.1415648) q[0];
sx q[0];
rz(0.4719032) q[0];
rz(0.42981237) q[1];
sx q[1];
rz(-1.2496354) q[1];
sx q[1];
rz(-0.93634161) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82114007) q[0];
sx q[0];
rz(-1.3409412) q[0];
sx q[0];
rz(-1.4810522) q[0];
rz(-0.78511946) q[2];
sx q[2];
rz(-1.8765419) q[2];
sx q[2];
rz(0.53346764) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1658926) q[1];
sx q[1];
rz(-0.70043889) q[1];
sx q[1];
rz(-2.9794934) q[1];
x q[2];
rz(2.0735047) q[3];
sx q[3];
rz(-1.5503251) q[3];
sx q[3];
rz(2.3308144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3669746) q[2];
sx q[2];
rz(-2.8145511) q[2];
sx q[2];
rz(0.42638391) q[2];
rz(1.9042227) q[3];
sx q[3];
rz(-0.62785134) q[3];
sx q[3];
rz(-3.1085076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24580978) q[0];
sx q[0];
rz(-1.3170467) q[0];
sx q[0];
rz(2.202503) q[0];
rz(0.89871961) q[1];
sx q[1];
rz(-0.4788613) q[1];
sx q[1];
rz(0.59392196) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31375162) q[0];
sx q[0];
rz(-1.8694287) q[0];
sx q[0];
rz(-1.9065501) q[0];
rz(-pi) q[1];
x q[1];
rz(0.057815876) q[2];
sx q[2];
rz(-0.88314344) q[2];
sx q[2];
rz(2.2857894) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2908823) q[1];
sx q[1];
rz(-1.1516654) q[1];
sx q[1];
rz(-0.23371975) q[1];
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
rz(2.5014191) q[2];
sx q[2];
rz(-0.72945014) q[2];
sx q[2];
rz(1.7017986) q[2];
rz(2.7539608) q[3];
sx q[3];
rz(-1.6250316) q[3];
sx q[3];
rz(-2.7534527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3751635) q[0];
sx q[0];
rz(-1.5901934) q[0];
sx q[0];
rz(-0.50278062) q[0];
rz(0.76820961) q[1];
sx q[1];
rz(-2.6380824) q[1];
sx q[1];
rz(2.3847413) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2198974) q[0];
sx q[0];
rz(-1.9031525) q[0];
sx q[0];
rz(2.7601348) q[0];
rz(0.012979522) q[2];
sx q[2];
rz(-2.0364967) q[2];
sx q[2];
rz(2.4109858) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.56843578) q[1];
sx q[1];
rz(-2.4385298) q[1];
sx q[1];
rz(-0.99848024) q[1];
rz(-pi) q[2];
rz(-0.73369153) q[3];
sx q[3];
rz(-1.1847704) q[3];
sx q[3];
rz(1.8611849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.42671529) q[2];
sx q[2];
rz(-1.2299512) q[2];
sx q[2];
rz(-1.4871917) q[2];
rz(-2.5590844) q[3];
sx q[3];
rz(-2.0472066) q[3];
sx q[3];
rz(2.5845161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29397598) q[0];
sx q[0];
rz(-1.0852381) q[0];
sx q[0];
rz(-0.75772444) q[0];
rz(1.853653) q[1];
sx q[1];
rz(-2.2133591) q[1];
sx q[1];
rz(-1.0505189) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30848635) q[0];
sx q[0];
rz(-1.3163438) q[0];
sx q[0];
rz(-1.893505) q[0];
x q[1];
rz(0.97425766) q[2];
sx q[2];
rz(-1.1968808) q[2];
sx q[2];
rz(-2.8982382) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1861021) q[1];
sx q[1];
rz(-1.0483861) q[1];
sx q[1];
rz(-2.2239457) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1220783) q[3];
sx q[3];
rz(-2.2973804) q[3];
sx q[3];
rz(-2.7725078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.22333764) q[2];
sx q[2];
rz(-0.35327062) q[2];
sx q[2];
rz(-2.5081432) q[2];
rz(-1.9472306) q[3];
sx q[3];
rz(-1.4712237) q[3];
sx q[3];
rz(0.71715322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
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
rz(-2.6514335) q[0];
sx q[0];
rz(-2.8884086) q[0];
rz(-1.6075915) q[1];
sx q[1];
rz(-1.4350767) q[1];
sx q[1];
rz(-1.4621428) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9896478) q[0];
sx q[0];
rz(-1.4914574) q[0];
sx q[0];
rz(2.9304855) q[0];
x q[1];
rz(-2.2199549) q[2];
sx q[2];
rz(-1.3772794) q[2];
sx q[2];
rz(2.9606539) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.51965442) q[1];
sx q[1];
rz(-2.7192273) q[1];
sx q[1];
rz(-0.42610355) q[1];
rz(-2.3092689) q[3];
sx q[3];
rz(-2.6905077) q[3];
sx q[3];
rz(-3.0576599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.2016466) q[2];
sx q[2];
rz(-0.74735171) q[2];
sx q[2];
rz(-2.3366826) q[2];
rz(1.9645875) q[3];
sx q[3];
rz(-1.0369119) q[3];
sx q[3];
rz(1.0837519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2729623) q[0];
sx q[0];
rz(-1.0694163) q[0];
sx q[0];
rz(-2.2139363) q[0];
rz(1.0246798) q[1];
sx q[1];
rz(-1.6352446) q[1];
sx q[1];
rz(2.129508) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1597848) q[0];
sx q[0];
rz(-0.73043981) q[0];
sx q[0];
rz(2.3083789) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.38883932) q[2];
sx q[2];
rz(-1.6992339) q[2];
sx q[2];
rz(-1.5156137) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.55793563) q[1];
sx q[1];
rz(-2.5968938) q[1];
sx q[1];
rz(-0.063636585) q[1];
rz(-pi) q[2];
rz(0.13109644) q[3];
sx q[3];
rz(-2.5790865) q[3];
sx q[3];
rz(1.8545811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6107789) q[2];
sx q[2];
rz(-1.6616219) q[2];
sx q[2];
rz(-0.42993316) q[2];
rz(2.1271465) q[3];
sx q[3];
rz(-0.40922624) q[3];
sx q[3];
rz(-2.608192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6124509) q[0];
sx q[0];
rz(-2.2021459) q[0];
sx q[0];
rz(-0.21417831) q[0];
rz(2.0902436) q[1];
sx q[1];
rz(-0.21251692) q[1];
sx q[1];
rz(2.8578551) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8026233) q[0];
sx q[0];
rz(-1.5366652) q[0];
sx q[0];
rz(-0.3011093) q[0];
x q[1];
rz(1.9010504) q[2];
sx q[2];
rz(-1.3888161) q[2];
sx q[2];
rz(-0.44588003) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6730496) q[1];
sx q[1];
rz(-1.0769516) q[1];
sx q[1];
rz(-1.3480575) q[1];
rz(-pi) q[2];
x q[2];
rz(0.65225668) q[3];
sx q[3];
rz(-1.0271003) q[3];
sx q[3];
rz(0.039747681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6909137) q[2];
sx q[2];
rz(-0.51270715) q[2];
sx q[2];
rz(1.696375) q[2];
rz(-1.5444267) q[3];
sx q[3];
rz(-1.701136) q[3];
sx q[3];
rz(-2.8022695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63672367) q[0];
sx q[0];
rz(-0.050014194) q[0];
sx q[0];
rz(-3.0723363) q[0];
rz(1.6537369) q[1];
sx q[1];
rz(-1.2830877) q[1];
sx q[1];
rz(1.5690631) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3587787) q[0];
sx q[0];
rz(-2.1920715) q[0];
sx q[0];
rz(-0.57399477) q[0];
rz(-pi) q[1];
rz(1.4885159) q[2];
sx q[2];
rz(-1.2459718) q[2];
sx q[2];
rz(2.5282853) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.442765) q[1];
sx q[1];
rz(-2.8932533) q[1];
sx q[1];
rz(0.34766867) q[1];
rz(-1.8435555) q[3];
sx q[3];
rz(-2.4747304) q[3];
sx q[3];
rz(-0.71344261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1853603) q[2];
sx q[2];
rz(-2.9113443) q[2];
sx q[2];
rz(0.13988477) q[2];
rz(-0.36758962) q[3];
sx q[3];
rz(-1.1871754) q[3];
sx q[3];
rz(2.1504413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1763828) q[0];
sx q[0];
rz(-2.7503224) q[0];
sx q[0];
rz(2.4998253) q[0];
rz(1.9104674) q[1];
sx q[1];
rz(-1.9893913) q[1];
sx q[1];
rz(-0.26783255) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5726686) q[0];
sx q[0];
rz(-2.2974282) q[0];
sx q[0];
rz(-2.5323244) q[0];
rz(1.7540625) q[2];
sx q[2];
rz(-1.6943309) q[2];
sx q[2];
rz(-2.0545517) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2512868) q[1];
sx q[1];
rz(-0.44154134) q[1];
sx q[1];
rz(-2.4114386) q[1];
rz(2.0249428) q[3];
sx q[3];
rz(-2.5513259) q[3];
sx q[3];
rz(0.86436194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3315167) q[2];
sx q[2];
rz(-1.6167567) q[2];
sx q[2];
rz(1.6646741) q[2];
rz(-0.26633513) q[3];
sx q[3];
rz(-0.24644066) q[3];
sx q[3];
rz(-0.5464856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1289566) q[0];
sx q[0];
rz(-0.92100443) q[0];
sx q[0];
rz(-2.2367649) q[0];
rz(0.77990445) q[1];
sx q[1];
rz(-2.654568) q[1];
sx q[1];
rz(1.754896) q[1];
rz(0.65200381) q[2];
sx q[2];
rz(-0.82521385) q[2];
sx q[2];
rz(3.0575976) q[2];
rz(-0.26364636) q[3];
sx q[3];
rz(-2.4468138) q[3];
sx q[3];
rz(3.1082603) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
