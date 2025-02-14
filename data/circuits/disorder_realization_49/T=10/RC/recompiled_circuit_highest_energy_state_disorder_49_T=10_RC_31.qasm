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
rz(0.3857412) q[0];
sx q[0];
rz(2.1585611) q[0];
sx q[0];
rz(9.9821363) q[0];
rz(1.214667) q[1];
sx q[1];
rz(-0.85433706) q[1];
sx q[1];
rz(2.1995423) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57948008) q[0];
sx q[0];
rz(-1.7796374) q[0];
sx q[0];
rz(-0.12128854) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4100907) q[2];
sx q[2];
rz(-2.3580551) q[2];
sx q[2];
rz(-2.6611626) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.97325402) q[1];
sx q[1];
rz(-1.1853519) q[1];
sx q[1];
rz(-1.5027352) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0090839) q[3];
sx q[3];
rz(-1.0367437) q[3];
sx q[3];
rz(-0.81919794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4152834) q[2];
sx q[2];
rz(-1.3741263) q[2];
sx q[2];
rz(1.7706002) q[2];
rz(-0.81909424) q[3];
sx q[3];
rz(-2.2947125) q[3];
sx q[3];
rz(0.11025652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46721989) q[0];
sx q[0];
rz(-1.233036) q[0];
sx q[0];
rz(-0.20587532) q[0];
rz(-1.8241833) q[1];
sx q[1];
rz(-1.2818047) q[1];
sx q[1];
rz(2.6726216) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7531573) q[0];
sx q[0];
rz(-2.4860365) q[0];
sx q[0];
rz(-3.1032733) q[0];
rz(-pi) q[1];
rz(-1.1577206) q[2];
sx q[2];
rz(-2.6402355) q[2];
sx q[2];
rz(2.219308) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4552942) q[1];
sx q[1];
rz(-2.1249692) q[1];
sx q[1];
rz(0.55880736) q[1];
rz(-pi) q[2];
rz(0.75403611) q[3];
sx q[3];
rz(-2.0181351) q[3];
sx q[3];
rz(3.0424117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.59490243) q[2];
sx q[2];
rz(-2.3089843) q[2];
sx q[2];
rz(-0.62758315) q[2];
rz(-2.0708496) q[3];
sx q[3];
rz(-2.6379733) q[3];
sx q[3];
rz(0.32881769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4048432) q[0];
sx q[0];
rz(-2.3312745) q[0];
sx q[0];
rz(2.3531083) q[0];
rz(2.5704747) q[1];
sx q[1];
rz(-1.8126789) q[1];
sx q[1];
rz(-1.4459389) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3841032) q[0];
sx q[0];
rz(-1.5840085) q[0];
sx q[0];
rz(1.8394952) q[0];
rz(-pi) q[1];
rz(0.064632434) q[2];
sx q[2];
rz(-2.6026313) q[2];
sx q[2];
rz(-0.85899437) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2647977) q[1];
sx q[1];
rz(-1.4709657) q[1];
sx q[1];
rz(-2.8431176) q[1];
rz(0.51579185) q[3];
sx q[3];
rz(-0.63871562) q[3];
sx q[3];
rz(-2.0206919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0995522) q[2];
sx q[2];
rz(-0.44241646) q[2];
sx q[2];
rz(-2.7624847) q[2];
rz(1.3742617) q[3];
sx q[3];
rz(-2.1702424) q[3];
sx q[3];
rz(2.1578535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9336201) q[0];
sx q[0];
rz(-2.3868028) q[0];
sx q[0];
rz(-2.2291613) q[0];
rz(-1.747793) q[1];
sx q[1];
rz(-2.6363966) q[1];
sx q[1];
rz(0.30232731) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1990406) q[0];
sx q[0];
rz(-1.4438757) q[0];
sx q[0];
rz(-1.5104483) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0197844) q[2];
sx q[2];
rz(-2.3752779) q[2];
sx q[2];
rz(-1.7036071) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6136377) q[1];
sx q[1];
rz(-0.21371811) q[1];
sx q[1];
rz(1.7760913) q[1];
rz(-pi) q[2];
rz(-2.264278) q[3];
sx q[3];
rz(-2.1319509) q[3];
sx q[3];
rz(-2.0356242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6572774) q[2];
sx q[2];
rz(-2.4399098) q[2];
sx q[2];
rz(-0.62937984) q[2];
rz(-1.7221919) q[3];
sx q[3];
rz(-1.281176) q[3];
sx q[3];
rz(2.4331376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0486384) q[0];
sx q[0];
rz(-2.104367) q[0];
sx q[0];
rz(-1.3193489) q[0];
rz(-2.0535779) q[1];
sx q[1];
rz(-1.8698147) q[1];
sx q[1];
rz(-1.6768657) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3528093) q[0];
sx q[0];
rz(-1.6029198) q[0];
sx q[0];
rz(0.87091586) q[0];
rz(2.508923) q[2];
sx q[2];
rz(-1.5969689) q[2];
sx q[2];
rz(2.2319792) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.32680997) q[1];
sx q[1];
rz(-1.9091354) q[1];
sx q[1];
rz(1.7717351) q[1];
rz(-pi) q[2];
rz(0.3220114) q[3];
sx q[3];
rz(-0.63068855) q[3];
sx q[3];
rz(1.2104152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4488843) q[2];
sx q[2];
rz(-3.0227737) q[2];
sx q[2];
rz(0.23692712) q[2];
rz(0.98053011) q[3];
sx q[3];
rz(-1.4292932) q[3];
sx q[3];
rz(2.197649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15810814) q[0];
sx q[0];
rz(-3.0321002) q[0];
sx q[0];
rz(-2.9964301) q[0];
rz(-1.6204087) q[1];
sx q[1];
rz(-2.3416134) q[1];
sx q[1];
rz(-3.1343585) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67381672) q[0];
sx q[0];
rz(-2.4002808) q[0];
sx q[0];
rz(-1.2178161) q[0];
x q[1];
rz(0.72609781) q[2];
sx q[2];
rz(-1.1396004) q[2];
sx q[2];
rz(0.83881718) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1177495) q[1];
sx q[1];
rz(-2.8850318) q[1];
sx q[1];
rz(1.9909977) q[1];
rz(-pi) q[2];
rz(1.5421914) q[3];
sx q[3];
rz(-1.4047457) q[3];
sx q[3];
rz(-0.19197026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.89473692) q[2];
sx q[2];
rz(-1.1892908) q[2];
sx q[2];
rz(2.5640633) q[2];
rz(-1.9891116) q[3];
sx q[3];
rz(-2.5566176) q[3];
sx q[3];
rz(1.729689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5685527) q[0];
sx q[0];
rz(-0.19565208) q[0];
sx q[0];
rz(-1.0618807) q[0];
rz(-1.7537687) q[1];
sx q[1];
rz(-1.9111218) q[1];
sx q[1];
rz(1.6434297) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1573054) q[0];
sx q[0];
rz(-1.1486736) q[0];
sx q[0];
rz(1.4253084) q[0];
x q[1];
rz(1.241356) q[2];
sx q[2];
rz(-0.95166517) q[2];
sx q[2];
rz(3.0513024) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5504341) q[1];
sx q[1];
rz(-0.60748749) q[1];
sx q[1];
rz(-1.1689069) q[1];
rz(-pi) q[2];
rz(0.56363799) q[3];
sx q[3];
rz(-2.5803714) q[3];
sx q[3];
rz(1.7983675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7428703) q[2];
sx q[2];
rz(-3.0169432) q[2];
sx q[2];
rz(-2.2275662) q[2];
rz(0.75383178) q[3];
sx q[3];
rz(-1.9099312) q[3];
sx q[3];
rz(2.6851173) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1330426) q[0];
sx q[0];
rz(-2.8164016) q[0];
sx q[0];
rz(3.131026) q[0];
rz(-2.1376624) q[1];
sx q[1];
rz(-1.7251451) q[1];
sx q[1];
rz(-3.1026057) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7332135) q[0];
sx q[0];
rz(-2.7762402) q[0];
sx q[0];
rz(1.062458) q[0];
rz(-pi) q[1];
rz(-2.8247897) q[2];
sx q[2];
rz(-1.5617019) q[2];
sx q[2];
rz(0.62745079) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7379563) q[1];
sx q[1];
rz(-1.9490598) q[1];
sx q[1];
rz(2.4191816) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.66308588) q[3];
sx q[3];
rz(-1.6414101) q[3];
sx q[3];
rz(-0.021573349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5799334) q[2];
sx q[2];
rz(-1.9114405) q[2];
sx q[2];
rz(3.0420493) q[2];
rz(1.8482515) q[3];
sx q[3];
rz(-1.2852531) q[3];
sx q[3];
rz(-0.37555638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.471591) q[0];
sx q[0];
rz(-1.3674068) q[0];
sx q[0];
rz(1.8294096) q[0];
rz(0.52681628) q[1];
sx q[1];
rz(-0.75378886) q[1];
sx q[1];
rz(-2.3048293) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.484641) q[0];
sx q[0];
rz(-1.8645727) q[0];
sx q[0];
rz(2.8696069) q[0];
rz(-pi) q[1];
rz(0.9545045) q[2];
sx q[2];
rz(-2.0319246) q[2];
sx q[2];
rz(-0.64660286) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2075314) q[1];
sx q[1];
rz(-2.4467) q[1];
sx q[1];
rz(2.7259105) q[1];
x q[2];
rz(-0.38133867) q[3];
sx q[3];
rz(-1.1313038) q[3];
sx q[3];
rz(2.5724996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3706563) q[2];
sx q[2];
rz(-0.61345658) q[2];
sx q[2];
rz(-1.9343617) q[2];
rz(1.3197445) q[3];
sx q[3];
rz(-1.5765669) q[3];
sx q[3];
rz(2.3453662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(0.094548263) q[0];
sx q[0];
rz(-0.47484174) q[0];
sx q[0];
rz(0.71665254) q[0];
rz(-1.563975) q[1];
sx q[1];
rz(-1.8378704) q[1];
sx q[1];
rz(0.61734739) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3302688) q[0];
sx q[0];
rz(-1.399938) q[0];
sx q[0];
rz(2.2854684) q[0];
rz(-pi) q[1];
rz(1.2812219) q[2];
sx q[2];
rz(-0.48557845) q[2];
sx q[2];
rz(-2.537279) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.58891728) q[1];
sx q[1];
rz(-2.008634) q[1];
sx q[1];
rz(-2.4312996) q[1];
rz(-pi) q[2];
rz(-0.93043296) q[3];
sx q[3];
rz(-1.798822) q[3];
sx q[3];
rz(-1.7785398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8120332) q[2];
sx q[2];
rz(-2.7541408) q[2];
sx q[2];
rz(2.6623902) q[2];
rz(1.6241578) q[3];
sx q[3];
rz(-1.4045709) q[3];
sx q[3];
rz(1.4994538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0064405) q[0];
sx q[0];
rz(-1.3545481) q[0];
sx q[0];
rz(-0.56726278) q[0];
rz(0.80815036) q[1];
sx q[1];
rz(-1.5222526) q[1];
sx q[1];
rz(1.6076988) q[1];
rz(-1.0508423) q[2];
sx q[2];
rz(-1.4563917) q[2];
sx q[2];
rz(0.2105486) q[2];
rz(0.80912374) q[3];
sx q[3];
rz(-1.3700784) q[3];
sx q[3];
rz(-2.6249052) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
