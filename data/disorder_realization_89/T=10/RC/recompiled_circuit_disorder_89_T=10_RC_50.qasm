OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.75582957) q[0];
sx q[0];
rz(-1.4094149) q[0];
sx q[0];
rz(-0.29456079) q[0];
rz(-2.8566868) q[1];
sx q[1];
rz(-2.6309738) q[1];
sx q[1];
rz(2.7198305) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.784214) q[0];
sx q[0];
rz(-1.5970018) q[0];
sx q[0];
rz(1.635301) q[0];
rz(-1.6099036) q[2];
sx q[2];
rz(-2.2570059) q[2];
sx q[2];
rz(1.4456911) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0114853) q[1];
sx q[1];
rz(-1.3890146) q[1];
sx q[1];
rz(-2.6891522) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5476417) q[3];
sx q[3];
rz(-1.6868601) q[3];
sx q[3];
rz(2.5719197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0191779) q[2];
sx q[2];
rz(-0.77029595) q[2];
sx q[2];
rz(2.1315234) q[2];
rz(1.4953556) q[3];
sx q[3];
rz(-1.3388747) q[3];
sx q[3];
rz(-8*pi/11) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0881969) q[0];
sx q[0];
rz(-1.2258376) q[0];
sx q[0];
rz(-2.2825867) q[0];
rz(2.7711218) q[1];
sx q[1];
rz(-1.5971239) q[1];
sx q[1];
rz(1.6765615) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1169491) q[0];
sx q[0];
rz(-1.878371) q[0];
sx q[0];
rz(1.3427539) q[0];
x q[1];
rz(0.53497603) q[2];
sx q[2];
rz(-1.2777395) q[2];
sx q[2];
rz(1.5154293) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4620004) q[1];
sx q[1];
rz(-0.781171) q[1];
sx q[1];
rz(-2.6189234) q[1];
rz(1.8584077) q[3];
sx q[3];
rz(-1.0606442) q[3];
sx q[3];
rz(0.19954296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7039965) q[2];
sx q[2];
rz(-1.2164755) q[2];
sx q[2];
rz(-0.55830467) q[2];
rz(-2.1022508) q[3];
sx q[3];
rz(-1.2149518) q[3];
sx q[3];
rz(-2.5487652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6168183) q[0];
sx q[0];
rz(-2.1891948) q[0];
sx q[0];
rz(-2.9779789) q[0];
rz(2.4257461) q[1];
sx q[1];
rz(-2.2952081) q[1];
sx q[1];
rz(1.3735501) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4931902) q[0];
sx q[0];
rz(-1.545854) q[0];
sx q[0];
rz(-0.52669749) q[0];
x q[1];
rz(2.2839196) q[2];
sx q[2];
rz(-1.3366941) q[2];
sx q[2];
rz(0.26299325) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9436685) q[1];
sx q[1];
rz(-1.0218044) q[1];
sx q[1];
rz(2.3180069) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3397066) q[3];
sx q[3];
rz(-0.31517866) q[3];
sx q[3];
rz(1.8568045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0107161) q[2];
sx q[2];
rz(-0.54543442) q[2];
sx q[2];
rz(-2.1219357) q[2];
rz(0.9807469) q[3];
sx q[3];
rz(-0.5766944) q[3];
sx q[3];
rz(-2.5286634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2402128) q[0];
sx q[0];
rz(-2.4432683) q[0];
sx q[0];
rz(2.3024978) q[0];
rz(-1.5377195) q[1];
sx q[1];
rz(-1.6042177) q[1];
sx q[1];
rz(2.531321) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82204098) q[0];
sx q[0];
rz(-0.10986957) q[0];
sx q[0];
rz(1.6739474) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7856366) q[2];
sx q[2];
rz(-1.9320556) q[2];
sx q[2];
rz(1.4796096) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0415045) q[1];
sx q[1];
rz(-1.3106292) q[1];
sx q[1];
rz(-2.2799904) q[1];
rz(-2.652076) q[3];
sx q[3];
rz(-2.094141) q[3];
sx q[3];
rz(2.948992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.443976) q[2];
sx q[2];
rz(-0.51091754) q[2];
sx q[2];
rz(-2.9042517) q[2];
rz(0.90732968) q[3];
sx q[3];
rz(-1.5932339) q[3];
sx q[3];
rz(1.2835519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5089371) q[0];
sx q[0];
rz(-0.11015686) q[0];
sx q[0];
rz(1.5326112) q[0];
rz(2.4081047) q[1];
sx q[1];
rz(-1.8746904) q[1];
sx q[1];
rz(-1.8283432) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26966306) q[0];
sx q[0];
rz(-0.73686826) q[0];
sx q[0];
rz(-2.065573) q[0];
rz(-pi) q[1];
rz(1.6275431) q[2];
sx q[2];
rz(-1.6065671) q[2];
sx q[2];
rz(-2.2030731) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.65391958) q[1];
sx q[1];
rz(-2.6012585) q[1];
sx q[1];
rz(0.78507702) q[1];
rz(-pi) q[2];
rz(-0.58532183) q[3];
sx q[3];
rz(-2.1748741) q[3];
sx q[3];
rz(0.6795336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0017172) q[2];
sx q[2];
rz(-1.2728609) q[2];
sx q[2];
rz(2.4772947) q[2];
rz(2.4364566) q[3];
sx q[3];
rz(-0.64283979) q[3];
sx q[3];
rz(0.10372182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.076684549) q[0];
sx q[0];
rz(-3.0397471) q[0];
sx q[0];
rz(-0.054811906) q[0];
rz(1.0143771) q[1];
sx q[1];
rz(-1.4784112) q[1];
sx q[1];
rz(0.48318133) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2830402) q[0];
sx q[0];
rz(-2.9498093) q[0];
sx q[0];
rz(1.900308) q[0];
rz(-pi) q[1];
rz(-1.6364355) q[2];
sx q[2];
rz(-0.46354957) q[2];
sx q[2];
rz(1.0077196) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.15258458) q[1];
sx q[1];
rz(-0.87859381) q[1];
sx q[1];
rz(-2.2425376) q[1];
rz(-pi) q[2];
rz(1.2504962) q[3];
sx q[3];
rz(-2.7926676) q[3];
sx q[3];
rz(2.2826209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.15423916) q[2];
sx q[2];
rz(-2.8581212) q[2];
sx q[2];
rz(3.0563291) q[2];
rz(-1.9412458) q[3];
sx q[3];
rz(-1.6253358) q[3];
sx q[3];
rz(2.9912662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36713704) q[0];
sx q[0];
rz(-0.13639233) q[0];
sx q[0];
rz(-0.95463395) q[0];
rz(-2.5700991) q[1];
sx q[1];
rz(-1.1936455) q[1];
sx q[1];
rz(-2.8894997) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.945767) q[0];
sx q[0];
rz(-0.44647631) q[0];
sx q[0];
rz(2.8004942) q[0];
rz(0.78105314) q[2];
sx q[2];
rz(-2.1456246) q[2];
sx q[2];
rz(2.1512254) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.39602173) q[1];
sx q[1];
rz(-0.77116291) q[1];
sx q[1];
rz(-1.039617) q[1];
rz(-pi) q[2];
rz(-2.0310568) q[3];
sx q[3];
rz(-1.9122951) q[3];
sx q[3];
rz(-0.5865435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7581042) q[2];
sx q[2];
rz(-2.1540116) q[2];
sx q[2];
rz(2.2371116) q[2];
rz(1.0036184) q[3];
sx q[3];
rz(-2.2763054) q[3];
sx q[3];
rz(-0.65902567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43338183) q[0];
sx q[0];
rz(-3.1383585) q[0];
sx q[0];
rz(-1.4138387) q[0];
rz(0.66043234) q[1];
sx q[1];
rz(-1.3884037) q[1];
sx q[1];
rz(-1.7699014) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46946457) q[0];
sx q[0];
rz(-0.95333316) q[0];
sx q[0];
rz(-3.0332964) q[0];
rz(0.22050627) q[2];
sx q[2];
rz(-0.77324152) q[2];
sx q[2];
rz(-2.3436433) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9904069) q[1];
sx q[1];
rz(-1.8670261) q[1];
sx q[1];
rz(-2.5469261) q[1];
rz(-pi) q[2];
rz(-1.4189818) q[3];
sx q[3];
rz(-1.9278798) q[3];
sx q[3];
rz(0.15541542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3147605) q[2];
sx q[2];
rz(-2.5173352) q[2];
sx q[2];
rz(-2.7436942) q[2];
rz(0.49063101) q[3];
sx q[3];
rz(-1.5450954) q[3];
sx q[3];
rz(-1.3354966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21815498) q[0];
sx q[0];
rz(-1.2175918) q[0];
sx q[0];
rz(2.9072705) q[0];
rz(1.6802457) q[1];
sx q[1];
rz(-0.82273465) q[1];
sx q[1];
rz(0.27059069) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9070248) q[0];
sx q[0];
rz(-1.7922635) q[0];
sx q[0];
rz(1.347876) q[0];
rz(1.9174689) q[2];
sx q[2];
rz(-1.7520095) q[2];
sx q[2];
rz(0.62435645) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2968263) q[1];
sx q[1];
rz(-1.4335732) q[1];
sx q[1];
rz(1.5674577) q[1];
x q[2];
rz(-1.3769334) q[3];
sx q[3];
rz(-0.48741515) q[3];
sx q[3];
rz(0.91538115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.33971912) q[2];
sx q[2];
rz(-0.084771307) q[2];
sx q[2];
rz(-1.6516997) q[2];
rz(2.4387032) q[3];
sx q[3];
rz(-1.1542164) q[3];
sx q[3];
rz(1.9809013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40437317) q[0];
sx q[0];
rz(-1.568856) q[0];
sx q[0];
rz(-2.9872966) q[0];
rz(2.1986296) q[1];
sx q[1];
rz(-1.1318726) q[1];
sx q[1];
rz(-2.399209) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61481793) q[0];
sx q[0];
rz(-0.97133884) q[0];
sx q[0];
rz(2.4012698) q[0];
rz(-2.7102091) q[2];
sx q[2];
rz(-2.1200392) q[2];
sx q[2];
rz(-2.0768349) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.91126373) q[1];
sx q[1];
rz(-1.2250369) q[1];
sx q[1];
rz(-1.5826349) q[1];
rz(0.14245716) q[3];
sx q[3];
rz(-1.5995306) q[3];
sx q[3];
rz(-1.8491668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8538889) q[2];
sx q[2];
rz(-2.0159857) q[2];
sx q[2];
rz(1.5819736) q[2];
rz(-0.21073267) q[3];
sx q[3];
rz(-0.78453523) q[3];
sx q[3];
rz(0.5334841) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29006526) q[0];
sx q[0];
rz(-1.0354488) q[0];
sx q[0];
rz(-2.5356472) q[0];
rz(1.3416946) q[1];
sx q[1];
rz(-1.1732027) q[1];
sx q[1];
rz(-1.8684594) q[1];
rz(-0.32158357) q[2];
sx q[2];
rz(-0.84761534) q[2];
sx q[2];
rz(-1.9615704) q[2];
rz(1.4349764) q[3];
sx q[3];
rz(-1.4321696) q[3];
sx q[3];
rz(-1.1925478) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];