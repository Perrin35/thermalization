OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.57920116) q[0];
sx q[0];
rz(5.5226749) q[0];
sx q[0];
rz(10.439846) q[0];
rz(1.9782344) q[1];
sx q[1];
rz(-0.42127633) q[1];
sx q[1];
rz(-1.9032698) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88841146) q[0];
sx q[0];
rz(-1.0602573) q[0];
sx q[0];
rz(-2.2967413) q[0];
rz(0.36739393) q[2];
sx q[2];
rz(-1.9695373) q[2];
sx q[2];
rz(-0.90786394) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.94436344) q[1];
sx q[1];
rz(-2.5355593) q[1];
sx q[1];
rz(-2.2418749) q[1];
rz(-pi) q[2];
rz(1.7389033) q[3];
sx q[3];
rz(-1.013275) q[3];
sx q[3];
rz(-1.3489189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6282661) q[2];
sx q[2];
rz(-1.3506177) q[2];
sx q[2];
rz(-0.70293054) q[2];
rz(-0.85567307) q[3];
sx q[3];
rz(-2.296505) q[3];
sx q[3];
rz(1.8446911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6913476) q[0];
sx q[0];
rz(-2.0424728) q[0];
sx q[0];
rz(0.74478373) q[0];
rz(1.567747) q[1];
sx q[1];
rz(-2.4796922) q[1];
sx q[1];
rz(-2.1994798) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73977913) q[0];
sx q[0];
rz(-1.3688068) q[0];
sx q[0];
rz(0.76633472) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6555384) q[2];
sx q[2];
rz(-2.3872445) q[2];
sx q[2];
rz(-0.95065439) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5444138) q[1];
sx q[1];
rz(-0.81477466) q[1];
sx q[1];
rz(0.18715231) q[1];
rz(-pi) q[2];
rz(-0.09047507) q[3];
sx q[3];
rz(-2.5535085) q[3];
sx q[3];
rz(-0.17234853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2524903) q[2];
sx q[2];
rz(-0.95688755) q[2];
sx q[2];
rz(-1.0955742) q[2];
rz(1.4787632) q[3];
sx q[3];
rz(-0.79083276) q[3];
sx q[3];
rz(1.3862632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11109322) q[0];
sx q[0];
rz(-1.4249304) q[0];
sx q[0];
rz(-1.7362562) q[0];
rz(1.7449215) q[1];
sx q[1];
rz(-1.7878572) q[1];
sx q[1];
rz(2.2415846) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8843061) q[0];
sx q[0];
rz(-2.2955756) q[0];
sx q[0];
rz(-2.0638564) q[0];
rz(-2.5609803) q[2];
sx q[2];
rz(-2.0502649) q[2];
sx q[2];
rz(1.1632077) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.72199539) q[1];
sx q[1];
rz(-1.743814) q[1];
sx q[1];
rz(0.74651511) q[1];
rz(-0.95696394) q[3];
sx q[3];
rz(-1.004289) q[3];
sx q[3];
rz(-1.5885799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6806543) q[2];
sx q[2];
rz(-2.9260981) q[2];
sx q[2];
rz(-0.94949618) q[2];
rz(2.5203868) q[3];
sx q[3];
rz(-2.216279) q[3];
sx q[3];
rz(1.1533823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42156521) q[0];
sx q[0];
rz(-1.8864487) q[0];
sx q[0];
rz(-0.048728745) q[0];
rz(0.85982927) q[1];
sx q[1];
rz(-0.33165926) q[1];
sx q[1];
rz(-2.8299832) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45845073) q[0];
sx q[0];
rz(-1.299017) q[0];
sx q[0];
rz(1.0807476) q[0];
rz(-pi) q[1];
rz(-2.7976102) q[2];
sx q[2];
rz(-2.0728353) q[2];
sx q[2];
rz(-3.0127522) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3032461) q[1];
sx q[1];
rz(-0.7184808) q[1];
sx q[1];
rz(-1.1619199) q[1];
rz(2.0120088) q[3];
sx q[3];
rz(-1.3008683) q[3];
sx q[3];
rz(0.91096151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9638046) q[2];
sx q[2];
rz(-0.21169855) q[2];
sx q[2];
rz(-1.6507899) q[2];
rz(-0.35495159) q[3];
sx q[3];
rz(-1.7345411) q[3];
sx q[3];
rz(1.2995592) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0519003) q[0];
sx q[0];
rz(-2.7841452) q[0];
sx q[0];
rz(1.8748913) q[0];
rz(-3.025324) q[1];
sx q[1];
rz(-0.99028844) q[1];
sx q[1];
rz(-1.6273392) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0851885) q[0];
sx q[0];
rz(-1.8833816) q[0];
sx q[0];
rz(0.87707918) q[0];
x q[1];
rz(-1.1852988) q[2];
sx q[2];
rz(-0.23699871) q[2];
sx q[2];
rz(-2.1182107) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0872495) q[1];
sx q[1];
rz(-1.5197481) q[1];
sx q[1];
rz(1.0264978) q[1];
rz(-pi) q[2];
rz(2.299304) q[3];
sx q[3];
rz(-1.5986048) q[3];
sx q[3];
rz(-2.9942346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2221471) q[2];
sx q[2];
rz(-2.6667892) q[2];
sx q[2];
rz(-1.1754645) q[2];
rz(0.90421024) q[3];
sx q[3];
rz(-1.892482) q[3];
sx q[3];
rz(-0.97833943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0745875) q[0];
sx q[0];
rz(-0.33127221) q[0];
sx q[0];
rz(-0.40147716) q[0];
rz(2.0145156) q[1];
sx q[1];
rz(-1.8311484) q[1];
sx q[1];
rz(2.3289767) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3298523) q[0];
sx q[0];
rz(-2.5074258) q[0];
sx q[0];
rz(-1.5638173) q[0];
rz(-pi) q[1];
rz(0.27490669) q[2];
sx q[2];
rz(-1.3454352) q[2];
sx q[2];
rz(1.7079086) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9795831) q[1];
sx q[1];
rz(-2.1623731) q[1];
sx q[1];
rz(1.3237557) q[1];
rz(-pi) q[2];
x q[2];
rz(0.62536247) q[3];
sx q[3];
rz(-2.2917622) q[3];
sx q[3];
rz(1.4537077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9408985) q[2];
sx q[2];
rz(-2.5914067) q[2];
sx q[2];
rz(1.5852488) q[2];
rz(0.19041348) q[3];
sx q[3];
rz(-0.75473458) q[3];
sx q[3];
rz(-1.2079027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3170526) q[0];
sx q[0];
rz(-0.38924488) q[0];
sx q[0];
rz(0.12271605) q[0];
rz(1.1931194) q[1];
sx q[1];
rz(-0.90843186) q[1];
sx q[1];
rz(0.76748031) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5106414) q[0];
sx q[0];
rz(-0.92483339) q[0];
sx q[0];
rz(-0.064152282) q[0];
x q[1];
rz(2.7597455) q[2];
sx q[2];
rz(-1.4918461) q[2];
sx q[2];
rz(-0.71982924) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2042522) q[1];
sx q[1];
rz(-2.2240727) q[1];
sx q[1];
rz(-3.0863347) q[1];
x q[2];
rz(0.010482739) q[3];
sx q[3];
rz(-0.40626486) q[3];
sx q[3];
rz(0.83028136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7439338) q[2];
sx q[2];
rz(-0.021641061) q[2];
sx q[2];
rz(-1.5896612) q[2];
rz(-1.6520366) q[3];
sx q[3];
rz(-1.7347615) q[3];
sx q[3];
rz(1.1154729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3234696) q[0];
sx q[0];
rz(-0.36190811) q[0];
sx q[0];
rz(-0.37539151) q[0];
rz(0.88343945) q[1];
sx q[1];
rz(-1.5366303) q[1];
sx q[1];
rz(-0.72867957) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8189296) q[0];
sx q[0];
rz(-2.0897341) q[0];
sx q[0];
rz(2.3170018) q[0];
rz(-pi) q[1];
x q[1];
rz(0.23587464) q[2];
sx q[2];
rz(-1.4067003) q[2];
sx q[2];
rz(2.1498888) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3672626) q[1];
sx q[1];
rz(-2.3891797) q[1];
sx q[1];
rz(1.03536) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1971134) q[3];
sx q[3];
rz(-0.91942838) q[3];
sx q[3];
rz(-2.8508972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4484619) q[2];
sx q[2];
rz(-1.5784266) q[2];
sx q[2];
rz(0.7473839) q[2];
rz(1.4061617) q[3];
sx q[3];
rz(-1.2179255) q[3];
sx q[3];
rz(-1.5555443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2324227) q[0];
sx q[0];
rz(-2.2012043) q[0];
sx q[0];
rz(0.7255834) q[0];
rz(1.8046509) q[1];
sx q[1];
rz(-1.2533816) q[1];
sx q[1];
rz(0.57428378) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0030356) q[0];
sx q[0];
rz(-1.4113562) q[0];
sx q[0];
rz(-0.11388679) q[0];
rz(-0.4094643) q[2];
sx q[2];
rz(-1.0396233) q[2];
sx q[2];
rz(-2.766618) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.060238801) q[1];
sx q[1];
rz(-1.4574086) q[1];
sx q[1];
rz(-2.3008623) q[1];
rz(-pi) q[2];
x q[2];
rz(0.78403715) q[3];
sx q[3];
rz(-1.136354) q[3];
sx q[3];
rz(1.8798352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1672704) q[2];
sx q[2];
rz(-1.7631301) q[2];
sx q[2];
rz(-0.66217011) q[2];
rz(2.3769489) q[3];
sx q[3];
rz(-1.60631) q[3];
sx q[3];
rz(-1.3242599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8793256) q[0];
sx q[0];
rz(-2.5536394) q[0];
sx q[0];
rz(-1.7437438) q[0];
rz(-1.0700048) q[1];
sx q[1];
rz(-2.013423) q[1];
sx q[1];
rz(-1.8096583) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71026826) q[0];
sx q[0];
rz(-1.6013751) q[0];
sx q[0];
rz(1.7215183) q[0];
x q[1];
rz(-1.9851763) q[2];
sx q[2];
rz(-1.7935026) q[2];
sx q[2];
rz(-3.0321995) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3146914) q[1];
sx q[1];
rz(-1.3413635) q[1];
sx q[1];
rz(1.8073842) q[1];
rz(-0.58140786) q[3];
sx q[3];
rz(-0.96145844) q[3];
sx q[3];
rz(1.4319624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3146882) q[2];
sx q[2];
rz(-2.0411699) q[2];
sx q[2];
rz(-0.17987128) q[2];
rz(1.3966857) q[3];
sx q[3];
rz(-1.7161918) q[3];
sx q[3];
rz(2.9250308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3708645) q[0];
sx q[0];
rz(-2.6612119) q[0];
sx q[0];
rz(-2.2330855) q[0];
rz(2.6188359) q[1];
sx q[1];
rz(-1.1154543) q[1];
sx q[1];
rz(1.9784068) q[1];
rz(-1.9942453) q[2];
sx q[2];
rz(-0.97931391) q[2];
sx q[2];
rz(-0.58936832) q[2];
rz(1.2046075) q[3];
sx q[3];
rz(-1.1259176) q[3];
sx q[3];
rz(2.068145) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
