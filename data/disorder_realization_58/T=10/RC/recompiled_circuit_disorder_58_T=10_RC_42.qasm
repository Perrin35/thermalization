OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9392202) q[0];
sx q[0];
rz(-0.4063172) q[0];
sx q[0];
rz(0.82011861) q[0];
rz(-0.36110538) q[1];
sx q[1];
rz(-2.5087924) q[1];
sx q[1];
rz(0.83067218) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2031189) q[0];
sx q[0];
rz(-1.1955368) q[0];
sx q[0];
rz(-1.9012326) q[0];
rz(2.4029998) q[2];
sx q[2];
rz(-1.2054218) q[2];
sx q[2];
rz(-1.5942758) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9484529) q[1];
sx q[1];
rz(-1.8857191) q[1];
sx q[1];
rz(0.83955168) q[1];
rz(-pi) q[2];
rz(-2.41483) q[3];
sx q[3];
rz(-0.57075497) q[3];
sx q[3];
rz(1.6970413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7044907) q[2];
sx q[2];
rz(-2.7097242) q[2];
sx q[2];
rz(-3.0214156) q[2];
rz(-1.9834571) q[3];
sx q[3];
rz(-1.3985671) q[3];
sx q[3];
rz(-0.91896287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.08081089) q[0];
sx q[0];
rz(-1.3588384) q[0];
sx q[0];
rz(-0.91180116) q[0];
rz(2.3520825) q[1];
sx q[1];
rz(-2.1531838) q[1];
sx q[1];
rz(2.8149014) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3347496) q[0];
sx q[0];
rz(-2.202889) q[0];
sx q[0];
rz(0.64081162) q[0];
rz(-1.2454883) q[2];
sx q[2];
rz(-2.6763958) q[2];
sx q[2];
rz(0.59782366) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.28194004) q[1];
sx q[1];
rz(-1.9743866) q[1];
sx q[1];
rz(2.2380026) q[1];
rz(-pi) q[2];
rz(3.0558415) q[3];
sx q[3];
rz(-2.2567856) q[3];
sx q[3];
rz(-3.0119197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.5102753) q[2];
sx q[2];
rz(-2.3687506) q[2];
sx q[2];
rz(-1.2871683) q[2];
rz(3.0316947) q[3];
sx q[3];
rz(-1.7307614) q[3];
sx q[3];
rz(-1.7597642) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.592955) q[0];
sx q[0];
rz(-2.4023963) q[0];
sx q[0];
rz(0.32989311) q[0];
rz(-2.864481) q[1];
sx q[1];
rz(-1.8246633) q[1];
sx q[1];
rz(-1.057391) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7085416) q[0];
sx q[0];
rz(-1.9110702) q[0];
sx q[0];
rz(-0.021854594) q[0];
rz(-pi) q[1];
rz(-1.3710652) q[2];
sx q[2];
rz(-1.9057416) q[2];
sx q[2];
rz(-1.0945601) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.33144618) q[1];
sx q[1];
rz(-2.3815037) q[1];
sx q[1];
rz(2.0600832) q[1];
rz(-pi) q[2];
rz(-2.7338545) q[3];
sx q[3];
rz(-0.45555112) q[3];
sx q[3];
rz(1.3726485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.7827591) q[2];
sx q[2];
rz(-1.1153355) q[2];
sx q[2];
rz(-1.7791629) q[2];
rz(0.6247012) q[3];
sx q[3];
rz(-1.0995068) q[3];
sx q[3];
rz(0.81956285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7841004) q[0];
sx q[0];
rz(-2.6153013) q[0];
sx q[0];
rz(-0.5350565) q[0];
rz(2.0013981) q[1];
sx q[1];
rz(-1.813872) q[1];
sx q[1];
rz(0.16539703) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3446942) q[0];
sx q[0];
rz(-1.7555153) q[0];
sx q[0];
rz(3.0152507) q[0];
rz(-0.16343127) q[2];
sx q[2];
rz(-1.7497517) q[2];
sx q[2];
rz(-0.93726678) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.71388984) q[1];
sx q[1];
rz(-2.1117003) q[1];
sx q[1];
rz(-1.9546024) q[1];
x q[2];
rz(1.0159675) q[3];
sx q[3];
rz(-1.8053683) q[3];
sx q[3];
rz(2.6382584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.39607221) q[2];
sx q[2];
rz(-1.8053651) q[2];
sx q[2];
rz(0.20544927) q[2];
rz(1.127634) q[3];
sx q[3];
rz(-1.1536359) q[3];
sx q[3];
rz(-1.2566465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4797392) q[0];
sx q[0];
rz(-2.2275708) q[0];
sx q[0];
rz(-2.7868295) q[0];
rz(1.154249) q[1];
sx q[1];
rz(-0.92461363) q[1];
sx q[1];
rz(0.23194557) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1922798) q[0];
sx q[0];
rz(-2.4711907) q[0];
sx q[0];
rz(-0.85286661) q[0];
x q[1];
rz(-0.16817981) q[2];
sx q[2];
rz(-0.7025223) q[2];
sx q[2];
rz(1.0934193) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.31084727) q[1];
sx q[1];
rz(-1.3535045) q[1];
sx q[1];
rz(-2.9160935) q[1];
rz(-1.5902953) q[3];
sx q[3];
rz(-1.9392506) q[3];
sx q[3];
rz(1.833545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.140124) q[2];
sx q[2];
rz(-1.3342369) q[2];
sx q[2];
rz(1.2403963) q[2];
rz(0.59605789) q[3];
sx q[3];
rz(-1.3052992) q[3];
sx q[3];
rz(1.8381455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0444788) q[0];
sx q[0];
rz(-3.0713186) q[0];
sx q[0];
rz(-0.20275673) q[0];
rz(0.98908201) q[1];
sx q[1];
rz(-1.443807) q[1];
sx q[1];
rz(2.1441377) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0925804) q[0];
sx q[0];
rz(-2.2708587) q[0];
sx q[0];
rz(-2.9655365) q[0];
rz(0.73156725) q[2];
sx q[2];
rz(-1.8533857) q[2];
sx q[2];
rz(0.76991316) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.62085405) q[1];
sx q[1];
rz(-1.8100909) q[1];
sx q[1];
rz(-2.005307) q[1];
rz(-pi) q[2];
x q[2];
rz(0.84647471) q[3];
sx q[3];
rz(-1.2292362) q[3];
sx q[3];
rz(-2.4404756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9138907) q[2];
sx q[2];
rz(-1.9753549) q[2];
sx q[2];
rz(0.4513936) q[2];
rz(0.40870062) q[3];
sx q[3];
rz(-1.5390076) q[3];
sx q[3];
rz(-1.9394978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8354427) q[0];
sx q[0];
rz(-2.7375484) q[0];
sx q[0];
rz(-0.62414449) q[0];
rz(1.5165326) q[1];
sx q[1];
rz(-0.25696483) q[1];
sx q[1];
rz(-0.61378941) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81543024) q[0];
sx q[0];
rz(-0.76061941) q[0];
sx q[0];
rz(1.8421696) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.21653793) q[2];
sx q[2];
rz(-2.8515186) q[2];
sx q[2];
rz(1.9921583) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.63749667) q[1];
sx q[1];
rz(-1.6747961) q[1];
sx q[1];
rz(1.6297479) q[1];
rz(-pi) q[2];
rz(-0.11717637) q[3];
sx q[3];
rz(-1.1083318) q[3];
sx q[3];
rz(-1.6346491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.43341407) q[2];
sx q[2];
rz(-2.2257979) q[2];
sx q[2];
rz(1.9667352) q[2];
rz(2.5332149) q[3];
sx q[3];
rz(-1.404168) q[3];
sx q[3];
rz(-1.7181989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2414918) q[0];
sx q[0];
rz(-1.8429723) q[0];
sx q[0];
rz(2.6532145) q[0];
rz(-1.5178559) q[1];
sx q[1];
rz(-1.7428215) q[1];
sx q[1];
rz(0.98446313) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2135293) q[0];
sx q[0];
rz(-2.5941879) q[0];
sx q[0];
rz(0.071896032) q[0];
rz(-pi) q[1];
rz(2.7266399) q[2];
sx q[2];
rz(-0.81595647) q[2];
sx q[2];
rz(0.50843898) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.24482803) q[1];
sx q[1];
rz(-2.7762189) q[1];
sx q[1];
rz(0.016896292) q[1];
rz(-pi) q[2];
rz(-1.8833313) q[3];
sx q[3];
rz(-2.0274037) q[3];
sx q[3];
rz(1.4708335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.70665923) q[2];
sx q[2];
rz(-1.443053) q[2];
sx q[2];
rz(-2.6064176) q[2];
rz(2.0914071) q[3];
sx q[3];
rz(-1.8015367) q[3];
sx q[3];
rz(-1.0092658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.500279) q[0];
sx q[0];
rz(-2.2135493) q[0];
sx q[0];
rz(2.4556659) q[0];
rz(-2.7507239) q[1];
sx q[1];
rz(-2.1978244) q[1];
sx q[1];
rz(2.2156782) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6322964) q[0];
sx q[0];
rz(-0.25521454) q[0];
sx q[0];
rz(-0.98310982) q[0];
x q[1];
rz(2.3472896) q[2];
sx q[2];
rz(-1.0216733) q[2];
sx q[2];
rz(1.6910451) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8402108) q[1];
sx q[1];
rz(-1.7173319) q[1];
sx q[1];
rz(-2.0135897) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4521452) q[3];
sx q[3];
rz(-1.3680397) q[3];
sx q[3];
rz(2.6076536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.19568504) q[2];
sx q[2];
rz(-2.6699799) q[2];
sx q[2];
rz(-0.53608981) q[2];
rz(2.7219971) q[3];
sx q[3];
rz(-2.5511238) q[3];
sx q[3];
rz(-1.6528116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0855899) q[0];
sx q[0];
rz(-0.35878006) q[0];
sx q[0];
rz(-2.7465903) q[0];
rz(1.6292054) q[1];
sx q[1];
rz(-1.2724266) q[1];
sx q[1];
rz(-1.013247) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17908827) q[0];
sx q[0];
rz(-0.27217406) q[0];
sx q[0];
rz(-2.1129235) q[0];
rz(-pi) q[1];
rz(-2.7878694) q[2];
sx q[2];
rz(-1.5295267) q[2];
sx q[2];
rz(1.8281787) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.80876795) q[1];
sx q[1];
rz(-1.1507532) q[1];
sx q[1];
rz(-1.4648449) q[1];
x q[2];
rz(0.89684422) q[3];
sx q[3];
rz(-2.2138811) q[3];
sx q[3];
rz(2.602592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.44832486) q[2];
sx q[2];
rz(-1.6222745) q[2];
sx q[2];
rz(-0.75941336) q[2];
rz(-1.7761207) q[3];
sx q[3];
rz(-2.0335734) q[3];
sx q[3];
rz(2.571648) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7286745) q[0];
sx q[0];
rz(-1.8771794) q[0];
sx q[0];
rz(-1.9369453) q[0];
rz(-0.57327523) q[1];
sx q[1];
rz(-1.234006) q[1];
sx q[1];
rz(-1.3201859) q[1];
rz(-1.4105994) q[2];
sx q[2];
rz(-0.44993958) q[2];
sx q[2];
rz(-1.1100563) q[2];
rz(-1.6462973) q[3];
sx q[3];
rz(-0.69926881) q[3];
sx q[3];
rz(0.26509501) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
