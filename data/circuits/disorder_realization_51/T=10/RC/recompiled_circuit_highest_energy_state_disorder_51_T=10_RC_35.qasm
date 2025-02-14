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
rz(0.32938862) q[0];
sx q[0];
rz(3.7731054) q[0];
sx q[0];
rz(9.4256529) q[0];
rz(2.5007091) q[1];
sx q[1];
rz(-2.1407318) q[1];
sx q[1];
rz(-0.34520087) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5986818) q[0];
sx q[0];
rz(-3.0474159) q[0];
sx q[0];
rz(1.3046632) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7183185) q[2];
sx q[2];
rz(-1.194724) q[2];
sx q[2];
rz(3.0126115) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7172437) q[1];
sx q[1];
rz(-2.1645438) q[1];
sx q[1];
rz(-2.6301066) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4963989) q[3];
sx q[3];
rz(-1.5870023) q[3];
sx q[3];
rz(-2.5287573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.22902809) q[2];
sx q[2];
rz(-1.3338858) q[2];
sx q[2];
rz(-2.933617) q[2];
rz(-2.8424756) q[3];
sx q[3];
rz(-0.57204539) q[3];
sx q[3];
rz(2.0007029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3308554) q[0];
sx q[0];
rz(-0.24092291) q[0];
sx q[0];
rz(-0.99288565) q[0];
rz(1.317124) q[1];
sx q[1];
rz(-2.8254852) q[1];
sx q[1];
rz(0.77450007) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74906384) q[0];
sx q[0];
rz(-1.3925902) q[0];
sx q[0];
rz(-3.130413) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.51807816) q[2];
sx q[2];
rz(-2.2642914) q[2];
sx q[2];
rz(0.77502807) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6933813) q[1];
sx q[1];
rz(-0.44084545) q[1];
sx q[1];
rz(-2.5930659) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9162972) q[3];
sx q[3];
rz(-1.1114128) q[3];
sx q[3];
rz(-0.53250203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.896686) q[2];
sx q[2];
rz(-1.8550355) q[2];
sx q[2];
rz(1.0150821) q[2];
rz(0.24909881) q[3];
sx q[3];
rz(-0.86779147) q[3];
sx q[3];
rz(-2.7740313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2772813) q[0];
sx q[0];
rz(-0.69121498) q[0];
sx q[0];
rz(-2.9840898) q[0];
rz(2.1306254) q[1];
sx q[1];
rz(-0.33282655) q[1];
sx q[1];
rz(-3.0027622) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1942753) q[0];
sx q[0];
rz(-1.2091769) q[0];
sx q[0];
rz(2.7951434) q[0];
rz(-pi) q[1];
rz(-2.6723318) q[2];
sx q[2];
rz(-0.77232328) q[2];
sx q[2];
rz(-2.0283716) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5032387) q[1];
sx q[1];
rz(-1.279083) q[1];
sx q[1];
rz(2.9364763) q[1];
rz(-pi) q[2];
rz(0.28335684) q[3];
sx q[3];
rz(-1.7390001) q[3];
sx q[3];
rz(-0.12780549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.46131721) q[2];
sx q[2];
rz(-2.2291144) q[2];
sx q[2];
rz(0.3581363) q[2];
rz(-0.67874587) q[3];
sx q[3];
rz(-0.94921422) q[3];
sx q[3];
rz(2.1024735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.045227483) q[0];
sx q[0];
rz(-2.1684833) q[0];
sx q[0];
rz(2.8367693) q[0];
rz(-0.40799704) q[1];
sx q[1];
rz(-1.7034737) q[1];
sx q[1];
rz(0.92794424) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1243195) q[0];
sx q[0];
rz(-1.7649179) q[0];
sx q[0];
rz(-0.13071901) q[0];
x q[1];
rz(0.60073845) q[2];
sx q[2];
rz(-2.1081807) q[2];
sx q[2];
rz(1.5058277) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2603307) q[1];
sx q[1];
rz(-1.7219543) q[1];
sx q[1];
rz(-2.4511454) q[1];
x q[2];
rz(-1.0545066) q[3];
sx q[3];
rz(-1.0060423) q[3];
sx q[3];
rz(1.3206583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1912332) q[2];
sx q[2];
rz(-0.57098907) q[2];
sx q[2];
rz(2.9534269) q[2];
rz(1.2763216) q[3];
sx q[3];
rz(-1.8264344) q[3];
sx q[3];
rz(1.7108542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4417878) q[0];
sx q[0];
rz(-0.34407523) q[0];
sx q[0];
rz(0.019388327) q[0];
rz(2.0806606) q[1];
sx q[1];
rz(-1.414199) q[1];
sx q[1];
rz(1.9020938) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0415619) q[0];
sx q[0];
rz(-1.2912371) q[0];
sx q[0];
rz(-0.97645219) q[0];
rz(-0.92186982) q[2];
sx q[2];
rz(-2.1989294) q[2];
sx q[2];
rz(0.79337304) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.46249786) q[1];
sx q[1];
rz(-2.0789008) q[1];
sx q[1];
rz(1.8930045) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.89770384) q[3];
sx q[3];
rz(-0.60808676) q[3];
sx q[3];
rz(-1.9185818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0197319) q[2];
sx q[2];
rz(-1.8283565) q[2];
sx q[2];
rz(1.8149553) q[2];
rz(2.895368) q[3];
sx q[3];
rz(-1.1028057) q[3];
sx q[3];
rz(-0.8518014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6061848) q[0];
sx q[0];
rz(-0.98537213) q[0];
sx q[0];
rz(-0.73915172) q[0];
rz(-2.7958561) q[1];
sx q[1];
rz(-1.812499) q[1];
sx q[1];
rz(0.87127042) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0381868) q[0];
sx q[0];
rz(-1.5015008) q[0];
sx q[0];
rz(2.359576) q[0];
rz(-pi) q[1];
rz(0.99314697) q[2];
sx q[2];
rz(-0.9075853) q[2];
sx q[2];
rz(3.1301067) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.024684357) q[1];
sx q[1];
rz(-1.4977807) q[1];
sx q[1];
rz(1.9597998) q[1];
x q[2];
rz(0.61136758) q[3];
sx q[3];
rz(-1.1725559) q[3];
sx q[3];
rz(-2.7643725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.31200108) q[2];
sx q[2];
rz(-2.5025554) q[2];
sx q[2];
rz(0.87257067) q[2];
rz(1.568659) q[3];
sx q[3];
rz(-0.60397732) q[3];
sx q[3];
rz(-0.93305552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0912112) q[0];
sx q[0];
rz(-2.8822883) q[0];
sx q[0];
rz(0.76164371) q[0];
rz(2.7161982) q[1];
sx q[1];
rz(-2.0014706) q[1];
sx q[1];
rz(-2.2147307) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.923553) q[0];
sx q[0];
rz(-1.3430376) q[0];
sx q[0];
rz(-0.24203288) q[0];
rz(-pi) q[1];
rz(-2.4285467) q[2];
sx q[2];
rz(-0.8536866) q[2];
sx q[2];
rz(2.6756659) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.95185223) q[1];
sx q[1];
rz(-1.9323255) q[1];
sx q[1];
rz(2.9041163) q[1];
rz(-pi) q[2];
rz(3.0934392) q[3];
sx q[3];
rz(-2.1934436) q[3];
sx q[3];
rz(-0.66648713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1130134) q[2];
sx q[2];
rz(-0.69602746) q[2];
sx q[2];
rz(-2.2704303) q[2];
rz(0.29843676) q[3];
sx q[3];
rz(-1.2454183) q[3];
sx q[3];
rz(1.3309853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4153862) q[0];
sx q[0];
rz(-2.6415249) q[0];
sx q[0];
rz(0.4183847) q[0];
rz(-0.2977953) q[1];
sx q[1];
rz(-1.9458385) q[1];
sx q[1];
rz(-3.0072838) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2373922) q[0];
sx q[0];
rz(-2.2456894) q[0];
sx q[0];
rz(-2.2740721) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9761132) q[2];
sx q[2];
rz(-1.7421054) q[2];
sx q[2];
rz(-1.5754981) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1301535) q[1];
sx q[1];
rz(-1.7053889) q[1];
sx q[1];
rz(2.9573739) q[1];
rz(-pi) q[2];
rz(-1.0487324) q[3];
sx q[3];
rz(-2.6122141) q[3];
sx q[3];
rz(0.10806187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.40992752) q[2];
sx q[2];
rz(-1.2890559) q[2];
sx q[2];
rz(-0.64201391) q[2];
rz(1.0632473) q[3];
sx q[3];
rz(-0.65170538) q[3];
sx q[3];
rz(-2.1003907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6006271) q[0];
sx q[0];
rz(-3.0525115) q[0];
sx q[0];
rz(-3.0294321) q[0];
rz(1.3345831) q[1];
sx q[1];
rz(-2.5490675) q[1];
sx q[1];
rz(-2.1122011) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5350069) q[0];
sx q[0];
rz(-2.2748053) q[0];
sx q[0];
rz(3.0993942) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.12597398) q[2];
sx q[2];
rz(-0.306189) q[2];
sx q[2];
rz(-0.15737113) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.411602) q[1];
sx q[1];
rz(-1.3407602) q[1];
sx q[1];
rz(3.0964628) q[1];
x q[2];
rz(1.4791895) q[3];
sx q[3];
rz(-2.0370968) q[3];
sx q[3];
rz(-2.8331851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.42334291) q[2];
sx q[2];
rz(-0.074904718) q[2];
sx q[2];
rz(-0.038012803) q[2];
rz(0.612261) q[3];
sx q[3];
rz(-2.2050048) q[3];
sx q[3];
rz(3.0003701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8925979) q[0];
sx q[0];
rz(-1.6706415) q[0];
sx q[0];
rz(2.7227962) q[0];
rz(-1.2576125) q[1];
sx q[1];
rz(-1.6323171) q[1];
sx q[1];
rz(-2.9990101) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52209848) q[0];
sx q[0];
rz(-1.5987248) q[0];
sx q[0];
rz(2.9832207) q[0];
rz(0.37303961) q[2];
sx q[2];
rz(-0.66290406) q[2];
sx q[2];
rz(2.474624) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0847454) q[1];
sx q[1];
rz(-1.0029625) q[1];
sx q[1];
rz(2.826087) q[1];
rz(2.7975818) q[3];
sx q[3];
rz(-1.1340525) q[3];
sx q[3];
rz(0.94319447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3451781) q[2];
sx q[2];
rz(-0.08064457) q[2];
sx q[2];
rz(0.26819116) q[2];
rz(-1.7372519) q[3];
sx q[3];
rz(-1.0920478) q[3];
sx q[3];
rz(-2.0973189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1254697) q[0];
sx q[0];
rz(-2.7653427) q[0];
sx q[0];
rz(0.34633037) q[0];
rz(-3.012433) q[1];
sx q[1];
rz(-1.2551413) q[1];
sx q[1];
rz(-1.7358949) q[1];
rz(1.2186981) q[2];
sx q[2];
rz(-1.1209956) q[2];
sx q[2];
rz(1.6605177) q[2];
rz(-3.075243) q[3];
sx q[3];
rz(-1.2766311) q[3];
sx q[3];
rz(1.0286812) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
