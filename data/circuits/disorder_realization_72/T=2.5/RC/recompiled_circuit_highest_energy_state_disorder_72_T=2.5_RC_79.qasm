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
rz(0.82062757) q[0];
sx q[0];
rz(-0.66250116) q[0];
sx q[0];
rz(-0.26354665) q[0];
rz(-2.0343556) q[1];
sx q[1];
rz(-1.2821481) q[1];
sx q[1];
rz(-0.67556226) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5476108) q[0];
sx q[0];
rz(-1.9063966) q[0];
sx q[0];
rz(-3.1347325) q[0];
rz(0.78111202) q[2];
sx q[2];
rz(-0.75133577) q[2];
sx q[2];
rz(-1.9066115) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.88982595) q[1];
sx q[1];
rz(-1.4205564) q[1];
sx q[1];
rz(1.9637693) q[1];
x q[2];
rz(1.0166024) q[3];
sx q[3];
rz(-2.5678291) q[3];
sx q[3];
rz(-2.2237908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.88966173) q[2];
sx q[2];
rz(-1.9193005) q[2];
sx q[2];
rz(-2.7737854) q[2];
rz(0.31428549) q[3];
sx q[3];
rz(-0.29044423) q[3];
sx q[3];
rz(-2.9634326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3322068) q[0];
sx q[0];
rz(-0.089501373) q[0];
sx q[0];
rz(1.7820763) q[0];
rz(1.2844194) q[1];
sx q[1];
rz(-1.9764683) q[1];
sx q[1];
rz(0.051609106) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8479659) q[0];
sx q[0];
rz(-0.32362263) q[0];
sx q[0];
rz(2.1725562) q[0];
rz(-0.20599761) q[2];
sx q[2];
rz(-1.9340746) q[2];
sx q[2];
rz(-1.8415123) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9549841) q[1];
sx q[1];
rz(-0.37395769) q[1];
sx q[1];
rz(-1.5723521) q[1];
rz(-pi) q[2];
rz(-2.6671131) q[3];
sx q[3];
rz(-1.9360376) q[3];
sx q[3];
rz(0.38680916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0525558) q[2];
sx q[2];
rz(-0.26036981) q[2];
sx q[2];
rz(2.5944749) q[2];
rz(1.268092) q[3];
sx q[3];
rz(-1.3601466) q[3];
sx q[3];
rz(-2.8030296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7608305) q[0];
sx q[0];
rz(-2.1050504) q[0];
sx q[0];
rz(-1.8007675) q[0];
rz(2.2876168) q[1];
sx q[1];
rz(-1.2869336) q[1];
sx q[1];
rz(-2.8391848) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9895565) q[0];
sx q[0];
rz(-1.2947417) q[0];
sx q[0];
rz(0.92854519) q[0];
rz(0.45291169) q[2];
sx q[2];
rz(-2.6476417) q[2];
sx q[2];
rz(-0.043722186) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1544309) q[1];
sx q[1];
rz(-2.2691474) q[1];
sx q[1];
rz(0.72583801) q[1];
rz(-pi) q[2];
rz(2.0624235) q[3];
sx q[3];
rz(-1.3280014) q[3];
sx q[3];
rz(-1.8427046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.93620187) q[2];
sx q[2];
rz(-1.256029) q[2];
sx q[2];
rz(-2.5993247) q[2];
rz(2.1451456) q[3];
sx q[3];
rz(-2.9139329) q[3];
sx q[3];
rz(-0.12797932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-2.9484613) q[0];
sx q[0];
rz(-1.6574991) q[0];
sx q[0];
rz(-1.2633854) q[0];
rz(2.7073233) q[1];
sx q[1];
rz(-2.3691005) q[1];
sx q[1];
rz(1.451452) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2094089) q[0];
sx q[0];
rz(-1.1157633) q[0];
sx q[0];
rz(-1.2463776) q[0];
rz(-2.473818) q[2];
sx q[2];
rz(-0.66168565) q[2];
sx q[2];
rz(-1.1061321) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5401104) q[1];
sx q[1];
rz(-0.14232351) q[1];
sx q[1];
rz(0.92592923) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2882388) q[3];
sx q[3];
rz(-2.382016) q[3];
sx q[3];
rz(-2.4977998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0749391) q[2];
sx q[2];
rz(-2.1982919) q[2];
sx q[2];
rz(2.2554876) q[2];
rz(-0.0063627176) q[3];
sx q[3];
rz(-1.3402091) q[3];
sx q[3];
rz(-1.1191012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5804382) q[0];
sx q[0];
rz(-0.4564603) q[0];
sx q[0];
rz(-1.2503257) q[0];
rz(-2.9087032) q[1];
sx q[1];
rz(-2.3857375) q[1];
sx q[1];
rz(-3.0511391) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6503432) q[0];
sx q[0];
rz(-1.3216265) q[0];
sx q[0];
rz(3.0155474) q[0];
rz(1.301342) q[2];
sx q[2];
rz(-1.1134992) q[2];
sx q[2];
rz(-0.85207613) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.0098086987) q[1];
sx q[1];
rz(-0.65612853) q[1];
sx q[1];
rz(-3.0387278) q[1];
rz(1.8188263) q[3];
sx q[3];
rz(-0.33735407) q[3];
sx q[3];
rz(-3.1243008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.49454409) q[2];
sx q[2];
rz(-2.9477951) q[2];
sx q[2];
rz(-1.9336644) q[2];
rz(0.9211933) q[3];
sx q[3];
rz(-2.2974206) q[3];
sx q[3];
rz(-2.6827961) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10814609) q[0];
sx q[0];
rz(-1.4220081) q[0];
sx q[0];
rz(-0.49149996) q[0];
rz(0.72832251) q[1];
sx q[1];
rz(-2.0305384) q[1];
sx q[1];
rz(-2.947015) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.174343) q[0];
sx q[0];
rz(-1.7922658) q[0];
sx q[0];
rz(-2.1465149) q[0];
rz(-pi) q[1];
rz(-2.3324703) q[2];
sx q[2];
rz(-2.1702907) q[2];
sx q[2];
rz(-0.15669151) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4763471) q[1];
sx q[1];
rz(-2.0788296) q[1];
sx q[1];
rz(-0.82583896) q[1];
rz(-0.65690144) q[3];
sx q[3];
rz(-1.3852775) q[3];
sx q[3];
rz(-0.82478648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3538094) q[2];
sx q[2];
rz(-2.1259191) q[2];
sx q[2];
rz(-0.0018399012) q[2];
rz(1.4932102) q[3];
sx q[3];
rz(-2.4108678) q[3];
sx q[3];
rz(-2.8572237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88347411) q[0];
sx q[0];
rz(-2.879877) q[0];
sx q[0];
rz(2.0765685) q[0];
rz(-0.26271543) q[1];
sx q[1];
rz(-0.5717259) q[1];
sx q[1];
rz(-0.46331847) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1048038) q[0];
sx q[0];
rz(-2.1533009) q[0];
sx q[0];
rz(-0.1881442) q[0];
rz(-2.5538374) q[2];
sx q[2];
rz(-0.81771933) q[2];
sx q[2];
rz(-2.9972434) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1541248) q[1];
sx q[1];
rz(-2.6374258) q[1];
sx q[1];
rz(2.8567863) q[1];
rz(-1.7780532) q[3];
sx q[3];
rz(-0.1175783) q[3];
sx q[3];
rz(-0.97299796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.00039438417) q[2];
sx q[2];
rz(-1.1606777) q[2];
sx q[2];
rz(-2.311643) q[2];
rz(-1.0478033) q[3];
sx q[3];
rz(-2.6959097) q[3];
sx q[3];
rz(0.1786264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8038427) q[0];
sx q[0];
rz(-2.1512845) q[0];
sx q[0];
rz(0.75482279) q[0];
rz(2.7524475) q[1];
sx q[1];
rz(-1.6837589) q[1];
sx q[1];
rz(0.39774242) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0254117) q[0];
sx q[0];
rz(-1.8278024) q[0];
sx q[0];
rz(-1.2776796) q[0];
rz(0.55126528) q[2];
sx q[2];
rz(-1.7659852) q[2];
sx q[2];
rz(-2.8697687) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.81521846) q[1];
sx q[1];
rz(-2.7644016) q[1];
sx q[1];
rz(-2.3227824) q[1];
rz(-pi) q[2];
x q[2];
rz(0.92614545) q[3];
sx q[3];
rz(-0.93346262) q[3];
sx q[3];
rz(1.6660575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8336746) q[2];
sx q[2];
rz(-0.29611823) q[2];
sx q[2];
rz(-0.31416848) q[2];
rz(1.2416154) q[3];
sx q[3];
rz(-1.5887458) q[3];
sx q[3];
rz(1.9023021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.817713) q[0];
sx q[0];
rz(-0.68809026) q[0];
sx q[0];
rz(2.896198) q[0];
rz(1.9112401) q[1];
sx q[1];
rz(-0.10086682) q[1];
sx q[1];
rz(2.6255677) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3115055) q[0];
sx q[0];
rz(-1.9507512) q[0];
sx q[0];
rz(1.7854693) q[0];
rz(-pi) q[1];
rz(2.4012309) q[2];
sx q[2];
rz(-0.90351331) q[2];
sx q[2];
rz(1.1297392) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.12052434) q[1];
sx q[1];
rz(-1.85441) q[1];
sx q[1];
rz(-2.7614158) q[1];
rz(-pi) q[2];
rz(1.6327758) q[3];
sx q[3];
rz(-1.4345285) q[3];
sx q[3];
rz(-0.34157066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5838991) q[2];
sx q[2];
rz(-1.2217174) q[2];
sx q[2];
rz(-2.058513) q[2];
rz(2.9205868) q[3];
sx q[3];
rz(-2.998304) q[3];
sx q[3];
rz(2.3451282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.021127064) q[0];
sx q[0];
rz(-0.080862232) q[0];
sx q[0];
rz(0.13588151) q[0];
rz(1.7120301) q[1];
sx q[1];
rz(-1.1212965) q[1];
sx q[1];
rz(0.28724614) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0025072) q[0];
sx q[0];
rz(-1.3268496) q[0];
sx q[0];
rz(0.19357817) q[0];
x q[1];
rz(-3.0570266) q[2];
sx q[2];
rz(-1.3340325) q[2];
sx q[2];
rz(2.6344476) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4212227) q[1];
sx q[1];
rz(-2.1541641) q[1];
sx q[1];
rz(-1.3066584) q[1];
x q[2];
rz(-2.8494037) q[3];
sx q[3];
rz(-1.0064408) q[3];
sx q[3];
rz(2.5404504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.91934943) q[2];
sx q[2];
rz(-1.9025849) q[2];
sx q[2];
rz(-0.83402056) q[2];
rz(-0.30759865) q[3];
sx q[3];
rz(-2.7643876) q[3];
sx q[3];
rz(2.8789177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1766227) q[0];
sx q[0];
rz(-1.3055834) q[0];
sx q[0];
rz(2.1156043) q[0];
rz(0.38656825) q[1];
sx q[1];
rz(-1.760512) q[1];
sx q[1];
rz(-0.61813933) q[1];
rz(-2.8564524) q[2];
sx q[2];
rz(-2.7013473) q[2];
sx q[2];
rz(-0.57821748) q[2];
rz(1.781842) q[3];
sx q[3];
rz(-1.4283709) q[3];
sx q[3];
rz(-1.4621468) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
