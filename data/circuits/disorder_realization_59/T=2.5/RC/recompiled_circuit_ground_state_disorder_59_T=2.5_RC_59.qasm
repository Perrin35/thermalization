OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0786809) q[0];
sx q[0];
rz(2.7432888) q[0];
sx q[0];
rz(9.7437133) q[0];
rz(2.6755264) q[1];
sx q[1];
rz(-1.4332486) q[1];
sx q[1];
rz(2.8679009) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6382054) q[0];
sx q[0];
rz(-2.2944258) q[0];
sx q[0];
rz(0.47576018) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.24342033) q[2];
sx q[2];
rz(-1.3491329) q[2];
sx q[2];
rz(-1.2482211) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4727386) q[1];
sx q[1];
rz(-1.3931119) q[1];
sx q[1];
rz(3.0230126) q[1];
x q[2];
rz(1.1496278) q[3];
sx q[3];
rz(-0.85794373) q[3];
sx q[3];
rz(1.1919392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8493001) q[2];
sx q[2];
rz(-0.92508525) q[2];
sx q[2];
rz(1.383847) q[2];
rz(-0.84550953) q[3];
sx q[3];
rz(-0.011761646) q[3];
sx q[3];
rz(-2.1498634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94075769) q[0];
sx q[0];
rz(-0.89460129) q[0];
sx q[0];
rz(1.0351329) q[0];
rz(-1.1677405) q[1];
sx q[1];
rz(-0.20226856) q[1];
sx q[1];
rz(-2.367173) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7073785) q[0];
sx q[0];
rz(-1.0814572) q[0];
sx q[0];
rz(2.6721832) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2317288) q[2];
sx q[2];
rz(-1.8468886) q[2];
sx q[2];
rz(-2.368449) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.32774156) q[1];
sx q[1];
rz(-2.6151513) q[1];
sx q[1];
rz(-1.3415643) q[1];
x q[2];
rz(-3.10765) q[3];
sx q[3];
rz(-1.5205624) q[3];
sx q[3];
rz(1.4242581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3356129) q[2];
sx q[2];
rz(-0.077294417) q[2];
sx q[2];
rz(0.94266164) q[2];
rz(3.0521657) q[3];
sx q[3];
rz(-0.663203) q[3];
sx q[3];
rz(2.8360143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5494004) q[0];
sx q[0];
rz(-0.58627272) q[0];
sx q[0];
rz(-0.19202448) q[0];
rz(-2.5281455) q[1];
sx q[1];
rz(-0.44812056) q[1];
sx q[1];
rz(0.014785756) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35280577) q[0];
sx q[0];
rz(-2.3067622) q[0];
sx q[0];
rz(2.1197181) q[0];
x q[1];
rz(-2.4528265) q[2];
sx q[2];
rz(-1.7851794) q[2];
sx q[2];
rz(2.2685879) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5303684) q[1];
sx q[1];
rz(-0.99677982) q[1];
sx q[1];
rz(-2.8635129) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9613033) q[3];
sx q[3];
rz(-3.0766682) q[3];
sx q[3];
rz(-2.8466259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.52745596) q[2];
sx q[2];
rz(-1.5946031) q[2];
sx q[2];
rz(3.1318437) q[2];
rz(2.5629432) q[3];
sx q[3];
rz(-2.0895683) q[3];
sx q[3];
rz(2.3600522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(3.0233199) q[0];
sx q[0];
rz(-2.3319722) q[0];
sx q[0];
rz(-0.6672346) q[0];
rz(1.0046593) q[1];
sx q[1];
rz(-0.15548448) q[1];
sx q[1];
rz(1.3803233) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8568958) q[0];
sx q[0];
rz(-2.0486437) q[0];
sx q[0];
rz(-0.010464593) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1491702) q[2];
sx q[2];
rz(-0.54900169) q[2];
sx q[2];
rz(-1.3379607) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8246824) q[1];
sx q[1];
rz(-1.0774634) q[1];
sx q[1];
rz(-1.6073336) q[1];
rz(-pi) q[2];
rz(2.0949239) q[3];
sx q[3];
rz(-1.7780684) q[3];
sx q[3];
rz(-0.85974271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.37463793) q[2];
sx q[2];
rz(-2.0071603) q[2];
sx q[2];
rz(-0.053442027) q[2];
rz(-2.5080569) q[3];
sx q[3];
rz(-0.41826785) q[3];
sx q[3];
rz(0.0088648349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87090129) q[0];
sx q[0];
rz(-0.57796657) q[0];
sx q[0];
rz(-2.0722678) q[0];
rz(-1.2879734) q[1];
sx q[1];
rz(-0.063871495) q[1];
sx q[1];
rz(3.0499444) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1081323) q[0];
sx q[0];
rz(-3.1164451) q[0];
sx q[0];
rz(2.0501686) q[0];
rz(2.1778278) q[2];
sx q[2];
rz(-0.53812611) q[2];
sx q[2];
rz(0.48689688) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7465072) q[1];
sx q[1];
rz(-1.2068492) q[1];
sx q[1];
rz(-2.0319315) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5318457) q[3];
sx q[3];
rz(-2.340303) q[3];
sx q[3];
rz(-1.0391446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2962467) q[2];
sx q[2];
rz(-2.3643934) q[2];
sx q[2];
rz(-0.11543154) q[2];
rz(2.3760065) q[3];
sx q[3];
rz(-1.6439532) q[3];
sx q[3];
rz(0.41281858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(2.5033443) q[0];
sx q[0];
rz(-0.64953506) q[0];
sx q[0];
rz(0.58393884) q[0];
rz(-3.0931547) q[1];
sx q[1];
rz(-0.2225288) q[1];
sx q[1];
rz(0.64360523) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2023425) q[0];
sx q[0];
rz(-1.5233375) q[0];
sx q[0];
rz(1.2623293) q[0];
x q[1];
rz(2.1435596) q[2];
sx q[2];
rz(-2.3946471) q[2];
sx q[2];
rz(-2.2446225) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.19932718) q[1];
sx q[1];
rz(-0.86546997) q[1];
sx q[1];
rz(-2.9678515) q[1];
x q[2];
rz(-2.3619283) q[3];
sx q[3];
rz(-1.7105967) q[3];
sx q[3];
rz(-1.9817022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.27759564) q[2];
sx q[2];
rz(-0.79825753) q[2];
sx q[2];
rz(-0.64211988) q[2];
rz(-3.0025205) q[3];
sx q[3];
rz(-0.13834794) q[3];
sx q[3];
rz(1.6819491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7794401) q[0];
sx q[0];
rz(-2.7208949) q[0];
sx q[0];
rz(0.62533373) q[0];
rz(-2.9026237) q[1];
sx q[1];
rz(-2.9049554) q[1];
sx q[1];
rz(1.3561358) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55955333) q[0];
sx q[0];
rz(-1.6252717) q[0];
sx q[0];
rz(-3.1097058) q[0];
x q[1];
rz(1.8618334) q[2];
sx q[2];
rz(-1.2087529) q[2];
sx q[2];
rz(-2.9933528) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6125638) q[1];
sx q[1];
rz(-0.41464889) q[1];
sx q[1];
rz(-1.0573483) q[1];
rz(-1.7271572) q[3];
sx q[3];
rz(-1.3605474) q[3];
sx q[3];
rz(2.4429136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.29505342) q[2];
sx q[2];
rz(-1.876538) q[2];
sx q[2];
rz(-0.96578252) q[2];
rz(0.24671181) q[3];
sx q[3];
rz(-1.8762981) q[3];
sx q[3];
rz(-1.3707967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5024289) q[0];
sx q[0];
rz(-0.40488365) q[0];
sx q[0];
rz(0.13614458) q[0];
rz(0.73774058) q[1];
sx q[1];
rz(-0.22519153) q[1];
sx q[1];
rz(-1.236261) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.284974) q[0];
sx q[0];
rz(-2.1559733) q[0];
sx q[0];
rz(-0.43584688) q[0];
x q[1];
rz(-1.3278373) q[2];
sx q[2];
rz(-2.4907673) q[2];
sx q[2];
rz(1.6215768) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7248316) q[1];
sx q[1];
rz(-0.5151075) q[1];
sx q[1];
rz(-2.5558042) q[1];
x q[2];
rz(1.1499722) q[3];
sx q[3];
rz(-0.69531402) q[3];
sx q[3];
rz(3.1041404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.69011921) q[2];
sx q[2];
rz(-1.5560919) q[2];
sx q[2];
rz(-2.1915009) q[2];
rz(-0.39725605) q[3];
sx q[3];
rz(-0.49644956) q[3];
sx q[3];
rz(-0.44601405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5600679) q[0];
sx q[0];
rz(-0.6811322) q[0];
sx q[0];
rz(0.11823046) q[0];
rz(-2.0589578) q[1];
sx q[1];
rz(-1.2018459) q[1];
sx q[1];
rz(-1.3523678) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0467606) q[0];
sx q[0];
rz(-2.7343395) q[0];
sx q[0];
rz(-2.880275) q[0];
rz(-1.6671343) q[2];
sx q[2];
rz(-2.1129449) q[2];
sx q[2];
rz(0.037329096) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0491421) q[1];
sx q[1];
rz(-2.051472) q[1];
sx q[1];
rz(2.9979826) q[1];
rz(-pi) q[2];
rz(1.0151789) q[3];
sx q[3];
rz(-2.5930282) q[3];
sx q[3];
rz(-0.11840222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3236986) q[2];
sx q[2];
rz(-2.4511621) q[2];
sx q[2];
rz(0.77198088) q[2];
rz(3.0585994) q[3];
sx q[3];
rz(-1.2677742) q[3];
sx q[3];
rz(0.52223372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8817187) q[0];
sx q[0];
rz(-2.3840955) q[0];
sx q[0];
rz(0.74257332) q[0];
rz(-1.2965797) q[1];
sx q[1];
rz(-2.3340338) q[1];
sx q[1];
rz(-0.47320941) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6542235) q[0];
sx q[0];
rz(-2.8608436) q[0];
sx q[0];
rz(1.8868179) q[0];
x q[1];
rz(2.8933018) q[2];
sx q[2];
rz(-1.8195229) q[2];
sx q[2];
rz(-1.7493713) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1291581) q[1];
sx q[1];
rz(-1.9170463) q[1];
sx q[1];
rz(-1.0770256) q[1];
rz(-pi) q[2];
x q[2];
rz(0.27373154) q[3];
sx q[3];
rz(-2.377284) q[3];
sx q[3];
rz(0.66018644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5903198) q[2];
sx q[2];
rz(-1.4164305) q[2];
sx q[2];
rz(1.9080706) q[2];
rz(2.7963855) q[3];
sx q[3];
rz(-0.39185169) q[3];
sx q[3];
rz(-2.3046618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45967669) q[0];
sx q[0];
rz(-1.7727333) q[0];
sx q[0];
rz(1.75417) q[0];
rz(0.058902901) q[1];
sx q[1];
rz(-1.7210996) q[1];
sx q[1];
rz(0.65239418) q[1];
rz(-2.7564213) q[2];
sx q[2];
rz(-1.9829911) q[2];
sx q[2];
rz(0.34172716) q[2];
rz(0.98109365) q[3];
sx q[3];
rz(-2.4193939) q[3];
sx q[3];
rz(-1.9388225) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
