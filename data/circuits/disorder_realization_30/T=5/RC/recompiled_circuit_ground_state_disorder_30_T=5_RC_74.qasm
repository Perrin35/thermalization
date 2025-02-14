OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8024017) q[0];
sx q[0];
rz(-1.985476) q[0];
sx q[0];
rz(-0.51751408) q[0];
rz(1.327688) q[1];
sx q[1];
rz(-2.2496536) q[1];
sx q[1];
rz(-2.5995624) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5919246) q[0];
sx q[0];
rz(-2.3506563) q[0];
sx q[0];
rz(0.24178453) q[0];
rz(0.26371358) q[2];
sx q[2];
rz(-0.78391916) q[2];
sx q[2];
rz(0.5413407) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.65329018) q[1];
sx q[1];
rz(-2.3494216) q[1];
sx q[1];
rz(1.8503055) q[1];
x q[2];
rz(-1.7446506) q[3];
sx q[3];
rz(-1.3411634) q[3];
sx q[3];
rz(-2.0876807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.0077670495) q[2];
sx q[2];
rz(-1.0453036) q[2];
sx q[2];
rz(-2.2110151) q[2];
rz(-3.0104356) q[3];
sx q[3];
rz(-2.7635062) q[3];
sx q[3];
rz(-1.650943) q[3];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59830484) q[0];
sx q[0];
rz(-0.85619339) q[0];
sx q[0];
rz(-1.244586) q[0];
rz(-2.2159131) q[1];
sx q[1];
rz(-1.2558179) q[1];
sx q[1];
rz(1.6474887) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3748532) q[0];
sx q[0];
rz(-1.245386) q[0];
sx q[0];
rz(-2.297657) q[0];
rz(-pi) q[1];
rz(-1.9060806) q[2];
sx q[2];
rz(-2.2163894) q[2];
sx q[2];
rz(-0.52926999) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0708965) q[1];
sx q[1];
rz(-0.58377171) q[1];
sx q[1];
rz(-2.4811005) q[1];
x q[2];
rz(3.1140987) q[3];
sx q[3];
rz(-0.3042694) q[3];
sx q[3];
rz(-2.9472443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1458448) q[2];
sx q[2];
rz(-2.6971942) q[2];
sx q[2];
rz(-0.90163976) q[2];
rz(1.7835279) q[3];
sx q[3];
rz(-1.8243022) q[3];
sx q[3];
rz(0.54734126) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77483344) q[0];
sx q[0];
rz(-1.1792553) q[0];
sx q[0];
rz(1.1460079) q[0];
rz(0.92578069) q[1];
sx q[1];
rz(-2.130276) q[1];
sx q[1];
rz(-2.944223) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5615047) q[0];
sx q[0];
rz(-0.60778032) q[0];
sx q[0];
rz(-1.32919) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8497265) q[2];
sx q[2];
rz(-1.2780683) q[2];
sx q[2];
rz(-1.9809993) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2132385) q[1];
sx q[1];
rz(-0.23121195) q[1];
sx q[1];
rz(-1.3998881) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3696543) q[3];
sx q[3];
rz(-2.1052868) q[3];
sx q[3];
rz(0.25132685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4027412) q[2];
sx q[2];
rz(-2.1612942) q[2];
sx q[2];
rz(2.951238) q[2];
rz(-0.061577408) q[3];
sx q[3];
rz(-0.96934167) q[3];
sx q[3];
rz(2.0891345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15678081) q[0];
sx q[0];
rz(-1.4182014) q[0];
sx q[0];
rz(0.62830997) q[0];
rz(1.620694) q[1];
sx q[1];
rz(-1.9338806) q[1];
sx q[1];
rz(0.77888387) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.395508) q[0];
sx q[0];
rz(-1.9451358) q[0];
sx q[0];
rz(0.94818481) q[0];
rz(-pi) q[1];
x q[1];
rz(0.36036307) q[2];
sx q[2];
rz(-1.833263) q[2];
sx q[2];
rz(-1.9356464) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.35793909) q[1];
sx q[1];
rz(-0.049555819) q[1];
sx q[1];
rz(2.8663584) q[1];
rz(-pi) q[2];
rz(-2.2728659) q[3];
sx q[3];
rz(-1.0208289) q[3];
sx q[3];
rz(1.3438674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4817619) q[2];
sx q[2];
rz(-1.0783106) q[2];
sx q[2];
rz(-3.1415494) q[2];
rz(-2.0569885) q[3];
sx q[3];
rz(-1.1559887) q[3];
sx q[3];
rz(0.46561766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67149177) q[0];
sx q[0];
rz(-1.4441613) q[0];
sx q[0];
rz(1.2985562) q[0];
rz(0.57688722) q[1];
sx q[1];
rz(-1.2737609) q[1];
sx q[1];
rz(-0.44688046) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80432207) q[0];
sx q[0];
rz(-1.7511586) q[0];
sx q[0];
rz(0.65482636) q[0];
x q[1];
rz(1.8700897) q[2];
sx q[2];
rz(-1.4993877) q[2];
sx q[2];
rz(1.7091319) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.19805867) q[1];
sx q[1];
rz(-1.7732278) q[1];
sx q[1];
rz(1.979503) q[1];
x q[2];
rz(-0.18549149) q[3];
sx q[3];
rz(-1.5122454) q[3];
sx q[3];
rz(1.0492988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.4345066) q[2];
sx q[2];
rz(-0.27927566) q[2];
sx q[2];
rz(2.9217829) q[2];
rz(1.6541727) q[3];
sx q[3];
rz(-1.4498962) q[3];
sx q[3];
rz(2.4421104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18911067) q[0];
sx q[0];
rz(-0.43287745) q[0];
sx q[0];
rz(2.2440946) q[0];
rz(-1.3709566) q[1];
sx q[1];
rz(-1.0400583) q[1];
sx q[1];
rz(-2.2460361) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5286551) q[0];
sx q[0];
rz(-0.85451689) q[0];
sx q[0];
rz(-1.3837293) q[0];
rz(0.95889556) q[2];
sx q[2];
rz(-1.0493426) q[2];
sx q[2];
rz(-2.796891) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0296391) q[1];
sx q[1];
rz(-1.9432507) q[1];
sx q[1];
rz(-1.2122494) q[1];
x q[2];
rz(1.0895133) q[3];
sx q[3];
rz(-1.1466422) q[3];
sx q[3];
rz(1.7103575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.36876496) q[2];
sx q[2];
rz(-1.1026829) q[2];
sx q[2];
rz(-0.15711288) q[2];
rz(-0.67240063) q[3];
sx q[3];
rz(-1.7417358) q[3];
sx q[3];
rz(3.0290643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2812578) q[0];
sx q[0];
rz(-2.4761138) q[0];
sx q[0];
rz(1.459664) q[0];
rz(-0.36059391) q[1];
sx q[1];
rz(-1.4692042) q[1];
sx q[1];
rz(1.2999387) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10541815) q[0];
sx q[0];
rz(-1.6929394) q[0];
sx q[0];
rz(2.7226177) q[0];
rz(-0.60827651) q[2];
sx q[2];
rz(-0.71675863) q[2];
sx q[2];
rz(0.83502095) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0850306) q[1];
sx q[1];
rz(-1.6675648) q[1];
sx q[1];
rz(-2.8438049) q[1];
rz(-pi) q[2];
rz(1.1645369) q[3];
sx q[3];
rz(-2.1463369) q[3];
sx q[3];
rz(-1.6353324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.931687) q[2];
sx q[2];
rz(-1.0301215) q[2];
sx q[2];
rz(-0.21200655) q[2];
rz(-2.0906406) q[3];
sx q[3];
rz(-1.7651599) q[3];
sx q[3];
rz(-0.088137805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5365005) q[0];
sx q[0];
rz(-2.3516646) q[0];
sx q[0];
rz(-2.9686046) q[0];
rz(1.1146924) q[1];
sx q[1];
rz(-1.0429017) q[1];
sx q[1];
rz(-0.10985049) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45013406) q[0];
sx q[0];
rz(-0.5930674) q[0];
sx q[0];
rz(0.35045464) q[0];
x q[1];
rz(-0.82288701) q[2];
sx q[2];
rz(-1.3776861) q[2];
sx q[2];
rz(0.88179526) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6348656) q[1];
sx q[1];
rz(-0.30931979) q[1];
sx q[1];
rz(0.97815467) q[1];
x q[2];
rz(0.19211003) q[3];
sx q[3];
rz(-2.5787836) q[3];
sx q[3];
rz(-2.5272688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4032119) q[2];
sx q[2];
rz(-1.2028368) q[2];
sx q[2];
rz(0.66486764) q[2];
rz(0.46857771) q[3];
sx q[3];
rz(-1.6178308) q[3];
sx q[3];
rz(0.81336462) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.453603) q[0];
sx q[0];
rz(-2.5316694) q[0];
sx q[0];
rz(0.74158057) q[0];
rz(-1.1874416) q[1];
sx q[1];
rz(-1.6148184) q[1];
sx q[1];
rz(-0.56914079) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3305107) q[0];
sx q[0];
rz(-1.2449322) q[0];
sx q[0];
rz(-1.261607) q[0];
x q[1];
rz(-0.14573614) q[2];
sx q[2];
rz(-0.45782858) q[2];
sx q[2];
rz(-0.65953461) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5829825) q[1];
sx q[1];
rz(-2.4931968) q[1];
sx q[1];
rz(-2.3984539) q[1];
rz(-0.69066234) q[3];
sx q[3];
rz(-0.96513018) q[3];
sx q[3];
rz(2.8298135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6471214) q[2];
sx q[2];
rz(-2.1529866) q[2];
sx q[2];
rz(-0.76623255) q[2];
rz(-1.1726441) q[3];
sx q[3];
rz(-1.6021043) q[3];
sx q[3];
rz(2.9197781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-2.2072993) q[0];
sx q[0];
rz(-0.3485637) q[0];
sx q[0];
rz(2.9497414) q[0];
rz(2.3244997) q[1];
sx q[1];
rz(-2.8849738) q[1];
sx q[1];
rz(0.70768913) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5930847) q[0];
sx q[0];
rz(-1.8157209) q[0];
sx q[0];
rz(-1.3913416) q[0];
rz(2.0470263) q[2];
sx q[2];
rz(-0.64733887) q[2];
sx q[2];
rz(-2.8425384) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.62218828) q[1];
sx q[1];
rz(-0.71336245) q[1];
sx q[1];
rz(-0.98825561) q[1];
rz(-0.79962279) q[3];
sx q[3];
rz(-0.76086894) q[3];
sx q[3];
rz(-1.7590673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2231458) q[2];
sx q[2];
rz(-0.94474363) q[2];
sx q[2];
rz(-0.72500149) q[2];
rz(-2.9265192) q[3];
sx q[3];
rz(-0.30078617) q[3];
sx q[3];
rz(2.422629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.219915) q[0];
sx q[0];
rz(-0.81659962) q[0];
sx q[0];
rz(0.30246977) q[0];
rz(-2.0401781) q[1];
sx q[1];
rz(-1.9093724) q[1];
sx q[1];
rz(-2.2473635) q[1];
rz(2.106059) q[2];
sx q[2];
rz(-1.4695833) q[2];
sx q[2];
rz(1.5493456) q[2];
rz(-0.054417944) q[3];
sx q[3];
rz(-2.6734753) q[3];
sx q[3];
rz(-2.5772167) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
