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
rz(1.6990868) q[0];
sx q[0];
rz(-1.8449755) q[0];
sx q[0];
rz(1.4427286) q[0];
rz(2.7891085) q[1];
sx q[1];
rz(-2.6152857) q[1];
sx q[1];
rz(-1.7736645) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57940021) q[0];
sx q[0];
rz(-2.693667) q[0];
sx q[0];
rz(2.5125487) q[0];
rz(-pi) q[1];
rz(-2.7858911) q[2];
sx q[2];
rz(-1.6061602) q[2];
sx q[2];
rz(0.48330467) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.38851967) q[1];
sx q[1];
rz(-2.6021829) q[1];
sx q[1];
rz(-2.8864157) q[1];
rz(-pi) q[2];
rz(2.519415) q[3];
sx q[3];
rz(-2.2832979) q[3];
sx q[3];
rz(-1.5714558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0679396) q[2];
sx q[2];
rz(-1.0509793) q[2];
sx q[2];
rz(-1.937872) q[2];
rz(0.97483557) q[3];
sx q[3];
rz(-0.9674415) q[3];
sx q[3];
rz(-1.8956641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3881653) q[0];
sx q[0];
rz(-1.0618671) q[0];
sx q[0];
rz(-0.88428307) q[0];
rz(0.70941225) q[1];
sx q[1];
rz(-1.8953036) q[1];
sx q[1];
rz(2.2394004) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1195148) q[0];
sx q[0];
rz(-1.365322) q[0];
sx q[0];
rz(-0.46066649) q[0];
x q[1];
rz(-2.1778267) q[2];
sx q[2];
rz(-1.6554518) q[2];
sx q[2];
rz(-0.36106685) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.480011) q[1];
sx q[1];
rz(-1.8152092) q[1];
sx q[1];
rz(2.8339425) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7092) q[3];
sx q[3];
rz(-1.4908199) q[3];
sx q[3];
rz(0.13025912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0238637) q[2];
sx q[2];
rz(-0.29568299) q[2];
sx q[2];
rz(1.4884865) q[2];
rz(-1.214341) q[3];
sx q[3];
rz(-1.4597273) q[3];
sx q[3];
rz(-1.2795718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.134267) q[0];
sx q[0];
rz(-1.3995582) q[0];
sx q[0];
rz(-0.16635995) q[0];
rz(0.69794377) q[1];
sx q[1];
rz(-2.3614387) q[1];
sx q[1];
rz(-1.4770329) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0255942) q[0];
sx q[0];
rz(-0.4420155) q[0];
sx q[0];
rz(1.5777977) q[0];
rz(-0.2715025) q[2];
sx q[2];
rz(-1.8418485) q[2];
sx q[2];
rz(2.0663911) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.046951357) q[1];
sx q[1];
rz(-1.1186386) q[1];
sx q[1];
rz(-0.96479123) q[1];
rz(-pi) q[2];
rz(0.03607492) q[3];
sx q[3];
rz(-0.49866184) q[3];
sx q[3];
rz(2.7447678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.6537689) q[2];
sx q[2];
rz(-1.1710125) q[2];
sx q[2];
rz(-0.6146532) q[2];
rz(1.6807618) q[3];
sx q[3];
rz(-1.5644282) q[3];
sx q[3];
rz(0.041725807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6093269) q[0];
sx q[0];
rz(-1.8121239) q[0];
sx q[0];
rz(-1.5186658) q[0];
rz(2.3329349) q[1];
sx q[1];
rz(-1.8452019) q[1];
sx q[1];
rz(3.0928968) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2653462) q[0];
sx q[0];
rz(-0.9101724) q[0];
sx q[0];
rz(-1.1424078) q[0];
rz(-pi) q[1];
rz(-0.28566177) q[2];
sx q[2];
rz(-2.7451519) q[2];
sx q[2];
rz(-0.91914058) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9295556) q[1];
sx q[1];
rz(-2.0494645) q[1];
sx q[1];
rz(-2.5655866) q[1];
rz(2.1118495) q[3];
sx q[3];
rz(-2.4623639) q[3];
sx q[3];
rz(1.7558869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0398756) q[2];
sx q[2];
rz(-2.3287435) q[2];
sx q[2];
rz(-1.8854878) q[2];
rz(-3.0806372) q[3];
sx q[3];
rz(-0.90164369) q[3];
sx q[3];
rz(-0.92756334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98546472) q[0];
sx q[0];
rz(-2.0138854) q[0];
sx q[0];
rz(-3.0319279) q[0];
rz(2.1127286) q[1];
sx q[1];
rz(-1.2867462) q[1];
sx q[1];
rz(1.5922155) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66707313) q[0];
sx q[0];
rz(-1.1231849) q[0];
sx q[0];
rz(-1.2569129) q[0];
x q[1];
rz(-2.778901) q[2];
sx q[2];
rz(-1.0981993) q[2];
sx q[2];
rz(-2.9031799) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2031496) q[1];
sx q[1];
rz(-1.1091653) q[1];
sx q[1];
rz(-2.0600956) q[1];
rz(1.6084592) q[3];
sx q[3];
rz(-1.6051014) q[3];
sx q[3];
rz(1.6543076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5900383) q[2];
sx q[2];
rz(-2.8503214) q[2];
sx q[2];
rz(0.48198286) q[2];
rz(2.7216952) q[3];
sx q[3];
rz(-1.0651383) q[3];
sx q[3];
rz(-0.70986748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14966203) q[0];
sx q[0];
rz(-0.88352942) q[0];
sx q[0];
rz(-0.71257198) q[0];
rz(-0.86992162) q[1];
sx q[1];
rz(-0.37418071) q[1];
sx q[1];
rz(3.1178927) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4159235) q[0];
sx q[0];
rz(-2.5953889) q[0];
sx q[0];
rz(-0.31809455) q[0];
rz(-2.4392088) q[2];
sx q[2];
rz(-2.0310035) q[2];
sx q[2];
rz(-0.61379877) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.99147532) q[1];
sx q[1];
rz(-1.5722701) q[1];
sx q[1];
rz(1.7365843) q[1];
rz(-pi) q[2];
rz(1.3605484) q[3];
sx q[3];
rz(-2.0422404) q[3];
sx q[3];
rz(1.8690848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.80453834) q[2];
sx q[2];
rz(-2.4961175) q[2];
sx q[2];
rz(-1.1402593) q[2];
rz(1.987847) q[3];
sx q[3];
rz(-1.9269582) q[3];
sx q[3];
rz(-1.8203189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.425728) q[0];
sx q[0];
rz(-1.7540997) q[0];
sx q[0];
rz(0.36780372) q[0];
rz(-1.4999207) q[1];
sx q[1];
rz(-0.92253128) q[1];
sx q[1];
rz(2.0466764) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2832868) q[0];
sx q[0];
rz(-1.0890766) q[0];
sx q[0];
rz(2.879309) q[0];
rz(-pi) q[1];
rz(-0.20873834) q[2];
sx q[2];
rz(-0.56025592) q[2];
sx q[2];
rz(1.7605361) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3907539) q[1];
sx q[1];
rz(-1.3350643) q[1];
sx q[1];
rz(1.5890934) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5209909) q[3];
sx q[3];
rz(-2.7375855) q[3];
sx q[3];
rz(1.9086407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0992004) q[2];
sx q[2];
rz(-1.3513869) q[2];
sx q[2];
rz(-3.0858827) q[2];
rz(-2.4646711) q[3];
sx q[3];
rz(-0.61546314) q[3];
sx q[3];
rz(-2.2630579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1894576) q[0];
sx q[0];
rz(-0.5624693) q[0];
sx q[0];
rz(-2.1761555) q[0];
rz(-1.2646487) q[1];
sx q[1];
rz(-2.3955884) q[1];
sx q[1];
rz(-2.8146578) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.801689) q[0];
sx q[0];
rz(-1.9205695) q[0];
sx q[0];
rz(3.0236493) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1189402) q[2];
sx q[2];
rz(-2.6955397) q[2];
sx q[2];
rz(1.2747916) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2911243) q[1];
sx q[1];
rz(-1.3936491) q[1];
sx q[1];
rz(-1.7825104) q[1];
rz(-3.0013004) q[3];
sx q[3];
rz(-2.6888118) q[3];
sx q[3];
rz(2.0425532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7586729) q[2];
sx q[2];
rz(-2.4372673) q[2];
sx q[2];
rz(2.3455589) q[2];
rz(0.087513611) q[3];
sx q[3];
rz(-0.60860601) q[3];
sx q[3];
rz(-2.2039738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4336808) q[0];
sx q[0];
rz(-1.7681363) q[0];
sx q[0];
rz(-0.32980907) q[0];
rz(0.073313035) q[1];
sx q[1];
rz(-0.70711702) q[1];
sx q[1];
rz(2.0475533) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1739376) q[0];
sx q[0];
rz(-0.47220818) q[0];
sx q[0];
rz(-3.0595991) q[0];
x q[1];
rz(-1.7297001) q[2];
sx q[2];
rz(-0.4022809) q[2];
sx q[2];
rz(-2.2251468) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0702522) q[1];
sx q[1];
rz(-2.1161181) q[1];
sx q[1];
rz(-1.7252183) q[1];
rz(0.039765914) q[3];
sx q[3];
rz(-1.13969) q[3];
sx q[3];
rz(0.72663166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0426992) q[2];
sx q[2];
rz(-0.77184474) q[2];
sx q[2];
rz(-1.137255) q[2];
rz(-0.62816652) q[3];
sx q[3];
rz(-0.38574949) q[3];
sx q[3];
rz(1.0073193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4758509) q[0];
sx q[0];
rz(-2.6306212) q[0];
sx q[0];
rz(-1.092859) q[0];
rz(-1.2650371) q[1];
sx q[1];
rz(-1.8362412) q[1];
sx q[1];
rz(-2.1078033) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2726934) q[0];
sx q[0];
rz(-1.0995691) q[0];
sx q[0];
rz(-2.1687105) q[0];
rz(-pi) q[1];
x q[1];
rz(0.7598147) q[2];
sx q[2];
rz(-1.7489479) q[2];
sx q[2];
rz(2.9896328) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7180192) q[1];
sx q[1];
rz(-1.1597654) q[1];
sx q[1];
rz(2.4573106) q[1];
rz(-pi) q[2];
rz(-1.4321253) q[3];
sx q[3];
rz(-1.6059173) q[3];
sx q[3];
rz(0.39681527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8826302) q[2];
sx q[2];
rz(-0.7236824) q[2];
sx q[2];
rz(-2.9222729) q[2];
rz(-2.3389881) q[3];
sx q[3];
rz(-1.31253) q[3];
sx q[3];
rz(2.820211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5103067) q[0];
sx q[0];
rz(-1.318537) q[0];
sx q[0];
rz(1.0986811) q[0];
rz(-2.9639099) q[1];
sx q[1];
rz(-2.1251353) q[1];
sx q[1];
rz(2.9716117) q[1];
rz(-1.4953407) q[2];
sx q[2];
rz(-0.70388489) q[2];
sx q[2];
rz(2.1415785) q[2];
rz(-1.7614638) q[3];
sx q[3];
rz(-1.4545961) q[3];
sx q[3];
rz(-2.0300099) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
