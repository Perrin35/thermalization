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
rz(-0.29210583) q[0];
sx q[0];
rz(-0.6120683) q[0];
sx q[0];
rz(-3.0906313) q[0];
rz(1.2904957) q[1];
sx q[1];
rz(-1.2761071) q[1];
sx q[1];
rz(2.2955503) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48117366) q[0];
sx q[0];
rz(-1.2221029) q[0];
sx q[0];
rz(1.3997188) q[0];
rz(-0.92535352) q[2];
sx q[2];
rz(-2.6820289) q[2];
sx q[2];
rz(-1.991635) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5436094) q[1];
sx q[1];
rz(-2.5777237) q[1];
sx q[1];
rz(-0.51582526) q[1];
x q[2];
rz(-1.5791164) q[3];
sx q[3];
rz(-1.6826311) q[3];
sx q[3];
rz(-0.40598265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7840665) q[2];
sx q[2];
rz(-0.94782031) q[2];
sx q[2];
rz(0.86304647) q[2];
rz(-2.0283902) q[3];
sx q[3];
rz(-0.92156711) q[3];
sx q[3];
rz(0.79871261) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5384101) q[0];
sx q[0];
rz(-0.69284678) q[0];
sx q[0];
rz(-2.8811654) q[0];
rz(-2.8159091) q[1];
sx q[1];
rz(-2.1574056) q[1];
sx q[1];
rz(-0.30776417) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69410283) q[0];
sx q[0];
rz(-1.8304951) q[0];
sx q[0];
rz(-2.3070241) q[0];
rz(-pi) q[1];
rz(-0.57450104) q[2];
sx q[2];
rz(-0.5443474) q[2];
sx q[2];
rz(-2.8050204) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.968281) q[1];
sx q[1];
rz(-1.3303583) q[1];
sx q[1];
rz(3.0467826) q[1];
rz(-pi) q[2];
rz(3.1347389) q[3];
sx q[3];
rz(-2.2061976) q[3];
sx q[3];
rz(1.8346661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8327568) q[2];
sx q[2];
rz(-0.86897659) q[2];
sx q[2];
rz(2.5858509) q[2];
rz(-0.51221383) q[3];
sx q[3];
rz(-2.5930976) q[3];
sx q[3];
rz(2.6223124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.25278768) q[0];
sx q[0];
rz(-0.35211173) q[0];
sx q[0];
rz(-2.3922701) q[0];
rz(-2.0454171) q[1];
sx q[1];
rz(-1.3744033) q[1];
sx q[1];
rz(-2.7574417) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6644728) q[0];
sx q[0];
rz(-2.4944803) q[0];
sx q[0];
rz(0.44575925) q[0];
rz(2.5403156) q[2];
sx q[2];
rz(-1.3601175) q[2];
sx q[2];
rz(2.7803068) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.13680116) q[1];
sx q[1];
rz(-0.30711949) q[1];
sx q[1];
rz(-1.7476487) q[1];
rz(-pi) q[2];
rz(2.8504653) q[3];
sx q[3];
rz(-0.90727057) q[3];
sx q[3];
rz(0.71602589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.38203794) q[2];
sx q[2];
rz(-1.3052772) q[2];
sx q[2];
rz(-1.6273512) q[2];
rz(2.4115327) q[3];
sx q[3];
rz(-0.74159139) q[3];
sx q[3];
rz(-0.39456427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45383129) q[0];
sx q[0];
rz(-1.7173959) q[0];
sx q[0];
rz(-0.060039595) q[0];
rz(2.2278348) q[1];
sx q[1];
rz(-2.2303631) q[1];
sx q[1];
rz(-0.93889108) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37591463) q[0];
sx q[0];
rz(-1.6041099) q[0];
sx q[0];
rz(3.0785962) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0591776) q[2];
sx q[2];
rz(-1.6019099) q[2];
sx q[2];
rz(1.7331725) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.36870391) q[1];
sx q[1];
rz(-1.0828754) q[1];
sx q[1];
rz(0.027603961) q[1];
rz(-0.77620971) q[3];
sx q[3];
rz(-0.84899711) q[3];
sx q[3];
rz(1.2152023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1303723) q[2];
sx q[2];
rz(-2.8488686) q[2];
sx q[2];
rz(0.65227738) q[2];
rz(1.8782714) q[3];
sx q[3];
rz(-1.5247366) q[3];
sx q[3];
rz(1.9704845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0703099) q[0];
sx q[0];
rz(-1.902453) q[0];
sx q[0];
rz(-2.5791445) q[0];
rz(1.7930454) q[1];
sx q[1];
rz(-0.28678539) q[1];
sx q[1];
rz(2.5602692) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5452257) q[0];
sx q[0];
rz(-2.2196182) q[0];
sx q[0];
rz(0.56745402) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6069218) q[2];
sx q[2];
rz(-2.0310406) q[2];
sx q[2];
rz(-2.4789916) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6299675) q[1];
sx q[1];
rz(-1.2164348) q[1];
sx q[1];
rz(2.0767861) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.082959) q[3];
sx q[3];
rz(-2.4982493) q[3];
sx q[3];
rz(1.7789121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.694146) q[2];
sx q[2];
rz(-1.4299102) q[2];
sx q[2];
rz(0.68667975) q[2];
rz(-2.741559) q[3];
sx q[3];
rz(-1.8828078) q[3];
sx q[3];
rz(0.011628477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.43468633) q[0];
sx q[0];
rz(-2.1676368) q[0];
sx q[0];
rz(2.6913225) q[0];
rz(-0.88092583) q[1];
sx q[1];
rz(-1.9073146) q[1];
sx q[1];
rz(1.1164104) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.085297708) q[0];
sx q[0];
rz(-1.3816815) q[0];
sx q[0];
rz(-2.3448639) q[0];
rz(-pi) q[1];
rz(-2.667465) q[2];
sx q[2];
rz(-0.38981405) q[2];
sx q[2];
rz(1.6451943) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7959545) q[1];
sx q[1];
rz(-1.5511723) q[1];
sx q[1];
rz(0.96052891) q[1];
rz(2.8003872) q[3];
sx q[3];
rz(-1.998105) q[3];
sx q[3];
rz(0.36580929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.66699666) q[2];
sx q[2];
rz(-2.907739) q[2];
sx q[2];
rz(2.4901566) q[2];
rz(0.78178072) q[3];
sx q[3];
rz(-1.4835446) q[3];
sx q[3];
rz(0.93530161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5930138) q[0];
sx q[0];
rz(-2.9055556) q[0];
sx q[0];
rz(0.4655984) q[0];
rz(1.1511401) q[1];
sx q[1];
rz(-1.2460037) q[1];
sx q[1];
rz(-1.3394855) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6205821) q[0];
sx q[0];
rz(-2.4537183) q[0];
sx q[0];
rz(2.4732711) q[0];
rz(-pi) q[1];
rz(2.2869682) q[2];
sx q[2];
rz(-3.0833088) q[2];
sx q[2];
rz(-0.56266498) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0420181) q[1];
sx q[1];
rz(-0.23222831) q[1];
sx q[1];
rz(-1.2378511) q[1];
x q[2];
rz(0.42134704) q[3];
sx q[3];
rz(-2.192492) q[3];
sx q[3];
rz(-0.58577116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0307978) q[2];
sx q[2];
rz(-1.75533) q[2];
sx q[2];
rz(3.1286855) q[2];
rz(-0.13478336) q[3];
sx q[3];
rz(-1.2603935) q[3];
sx q[3];
rz(2.9002262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(1.6808692) q[0];
sx q[0];
rz(-0.091982059) q[0];
sx q[0];
rz(1.1750093) q[0];
rz(2.3911632) q[1];
sx q[1];
rz(-1.78396) q[1];
sx q[1];
rz(0.39355412) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4884696) q[0];
sx q[0];
rz(-2.5284889) q[0];
sx q[0];
rz(-1.8777008) q[0];
rz(2.9904994) q[2];
sx q[2];
rz(-2.4799621) q[2];
sx q[2];
rz(2.2349295) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4992884) q[1];
sx q[1];
rz(-2.6975216) q[1];
sx q[1];
rz(-1.8288235) q[1];
rz(0.39799798) q[3];
sx q[3];
rz(-1.0915712) q[3];
sx q[3];
rz(2.3204539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5081818) q[2];
sx q[2];
rz(-1.1339302) q[2];
sx q[2];
rz(0.22937648) q[2];
rz(-3.0502099) q[3];
sx q[3];
rz(-1.6978076) q[3];
sx q[3];
rz(2.5858333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9230187) q[0];
sx q[0];
rz(-0.78892437) q[0];
sx q[0];
rz(2.5564585) q[0];
rz(1.2773889) q[1];
sx q[1];
rz(-1.9265415) q[1];
sx q[1];
rz(-2.9871984) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.383787) q[0];
sx q[0];
rz(-0.4497779) q[0];
sx q[0];
rz(-0.25346724) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2227258) q[2];
sx q[2];
rz(-1.1670975) q[2];
sx q[2];
rz(-2.8793546) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8081863) q[1];
sx q[1];
rz(-1.8318614) q[1];
sx q[1];
rz(0.33269791) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2326689) q[3];
sx q[3];
rz(-2.6498389) q[3];
sx q[3];
rz(-1.7724747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1207017) q[2];
sx q[2];
rz(-2.2521844) q[2];
sx q[2];
rz(0.54023877) q[2];
rz(0.39515105) q[3];
sx q[3];
rz(-2.4873147) q[3];
sx q[3];
rz(-1.7925709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9304792) q[0];
sx q[0];
rz(-1.2808639) q[0];
sx q[0];
rz(-1.5323918) q[0];
rz(2.2156175) q[1];
sx q[1];
rz(-1.4264868) q[1];
sx q[1];
rz(1.4769953) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84718448) q[0];
sx q[0];
rz(-1.5715181) q[0];
sx q[0];
rz(2.4154354) q[0];
rz(-pi) q[1];
rz(-1.9473929) q[2];
sx q[2];
rz(-1.8459326) q[2];
sx q[2];
rz(-1.1801495) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.033869318) q[1];
sx q[1];
rz(-2.1185912) q[1];
sx q[1];
rz(-2.7366927) q[1];
rz(-pi) q[2];
rz(-2.1595803) q[3];
sx q[3];
rz(-1.1945621) q[3];
sx q[3];
rz(-2.328095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.28107873) q[2];
sx q[2];
rz(-1.8033359) q[2];
sx q[2];
rz(-2.4930084) q[2];
rz(-0.72566882) q[3];
sx q[3];
rz(-0.23894335) q[3];
sx q[3];
rz(1.749595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3758748) q[0];
sx q[0];
rz(-2.1633042) q[0];
sx q[0];
rz(-2.3194585) q[0];
rz(-2.0769465) q[1];
sx q[1];
rz(-1.6864265) q[1];
sx q[1];
rz(2.0815157) q[1];
rz(-0.2343455) q[2];
sx q[2];
rz(-2.4916983) q[2];
sx q[2];
rz(-1.1532952) q[2];
rz(1.0242994) q[3];
sx q[3];
rz(-2.306675) q[3];
sx q[3];
rz(-2.6425895) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
