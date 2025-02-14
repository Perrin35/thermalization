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
rz(1.3069557) q[0];
sx q[0];
rz(-2.262158) q[0];
sx q[0];
rz(0.49361324) q[0];
rz(-1.2797132) q[1];
sx q[1];
rz(-2.5149829) q[1];
sx q[1];
rz(1.6630747) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2459697) q[0];
sx q[0];
rz(-0.9446656) q[0];
sx q[0];
rz(0.62762053) q[0];
x q[1];
rz(0.25716146) q[2];
sx q[2];
rz(-2.2182121) q[2];
sx q[2];
rz(1.1360628) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.263674) q[1];
sx q[1];
rz(-1.3527737) q[1];
sx q[1];
rz(0.71123567) q[1];
x q[2];
rz(-1.0888723) q[3];
sx q[3];
rz(-2.4471209) q[3];
sx q[3];
rz(-2.0884909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0066068) q[2];
sx q[2];
rz(-2.0015494) q[2];
sx q[2];
rz(2.0331649) q[2];
rz(-1.3125575) q[3];
sx q[3];
rz(-2.8969942) q[3];
sx q[3];
rz(0.31497064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38158622) q[0];
sx q[0];
rz(-0.8257603) q[0];
sx q[0];
rz(-0.015856892) q[0];
rz(-2.1511757) q[1];
sx q[1];
rz(-2.3742193) q[1];
sx q[1];
rz(-2.2784065) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6619157) q[0];
sx q[0];
rz(-0.74115314) q[0];
sx q[0];
rz(0.14314289) q[0];
rz(2.3658889) q[2];
sx q[2];
rz(-2.687907) q[2];
sx q[2];
rz(2.016748) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.51746115) q[1];
sx q[1];
rz(-1.7449813) q[1];
sx q[1];
rz(0.49887533) q[1];
rz(-pi) q[2];
rz(0.84896481) q[3];
sx q[3];
rz(-1.6833653) q[3];
sx q[3];
rz(1.6602885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4414759) q[2];
sx q[2];
rz(-0.2773383) q[2];
sx q[2];
rz(0.46879834) q[2];
rz(0.68764728) q[3];
sx q[3];
rz(-1.629849) q[3];
sx q[3];
rz(-2.8414753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.0013292) q[0];
sx q[0];
rz(-0.92324531) q[0];
sx q[0];
rz(-0.73626751) q[0];
rz(0.34321347) q[1];
sx q[1];
rz(-2.2540269) q[1];
sx q[1];
rz(0.18377486) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7385173) q[0];
sx q[0];
rz(-2.3949497) q[0];
sx q[0];
rz(3.0439506) q[0];
x q[1];
rz(1.0509346) q[2];
sx q[2];
rz(-0.94280548) q[2];
sx q[2];
rz(1.7597212) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6334875) q[1];
sx q[1];
rz(-2.268666) q[1];
sx q[1];
rz(-1.1771604) q[1];
rz(-pi) q[2];
rz(-2.1101489) q[3];
sx q[3];
rz(-1.8832409) q[3];
sx q[3];
rz(-0.38020089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3595714) q[2];
sx q[2];
rz(-2.6247793) q[2];
sx q[2];
rz(-0.51592958) q[2];
rz(2.0333717) q[3];
sx q[3];
rz(-0.55086946) q[3];
sx q[3];
rz(1.1512604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(2.9786523) q[0];
sx q[0];
rz(-1.9300224) q[0];
sx q[0];
rz(-0.51280713) q[0];
rz(2.6817952) q[1];
sx q[1];
rz(-1.5922981) q[1];
sx q[1];
rz(-2.3777681) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5483185) q[0];
sx q[0];
rz(-1.700907) q[0];
sx q[0];
rz(-1.4553372) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.41641367) q[2];
sx q[2];
rz(-1.8685307) q[2];
sx q[2];
rz(-2.1418051) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0415548) q[1];
sx q[1];
rz(-0.44719346) q[1];
sx q[1];
rz(-2.2639422) q[1];
rz(-pi) q[2];
rz(-2.7598792) q[3];
sx q[3];
rz(-0.36892051) q[3];
sx q[3];
rz(-2.5526508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1226471) q[2];
sx q[2];
rz(-2.985869) q[2];
sx q[2];
rz(-0.32749185) q[2];
rz(0.28199768) q[3];
sx q[3];
rz(-1.7433386) q[3];
sx q[3];
rz(1.5544844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9649488) q[0];
sx q[0];
rz(-2.7556941) q[0];
sx q[0];
rz(-0.97736812) q[0];
rz(0.36879677) q[1];
sx q[1];
rz(-2.2803523) q[1];
sx q[1];
rz(-1.7288953) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0127129) q[0];
sx q[0];
rz(-1.0815718) q[0];
sx q[0];
rz(0.21531824) q[0];
rz(-pi) q[1];
rz(2.0875889) q[2];
sx q[2];
rz(-1.5587423) q[2];
sx q[2];
rz(0.64753676) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1407847) q[1];
sx q[1];
rz(-0.28433263) q[1];
sx q[1];
rz(1.4592912) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8983973) q[3];
sx q[3];
rz(-2.1246111) q[3];
sx q[3];
rz(-0.49529759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.29350027) q[2];
sx q[2];
rz(-2.3257181) q[2];
sx q[2];
rz(2.9917742) q[2];
rz(-0.39854974) q[3];
sx q[3];
rz(-1.6179251) q[3];
sx q[3];
rz(-2.985305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77001101) q[0];
sx q[0];
rz(-0.12938975) q[0];
sx q[0];
rz(2.5883801) q[0];
rz(2.5999056) q[1];
sx q[1];
rz(-1.0463511) q[1];
sx q[1];
rz(0.4479301) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89066896) q[0];
sx q[0];
rz(-0.4145997) q[0];
sx q[0];
rz(2.7756734) q[0];
rz(-0.42945736) q[2];
sx q[2];
rz(-1.1633368) q[2];
sx q[2];
rz(-2.7384659) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6363931) q[1];
sx q[1];
rz(-1.4318887) q[1];
sx q[1];
rz(-1.1769151) q[1];
x q[2];
rz(2.9782441) q[3];
sx q[3];
rz(-0.33698296) q[3];
sx q[3];
rz(-0.28675045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.64122671) q[2];
sx q[2];
rz(-1.9772823) q[2];
sx q[2];
rz(0.79916239) q[2];
rz(2.7568119) q[3];
sx q[3];
rz(-2.81541) q[3];
sx q[3];
rz(1.0001812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7810829) q[0];
sx q[0];
rz(-2.0979083) q[0];
sx q[0];
rz(-0.59180301) q[0];
rz(2.7599755) q[1];
sx q[1];
rz(-2.7615669) q[1];
sx q[1];
rz(0.15242481) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8417609) q[0];
sx q[0];
rz(-2.526847) q[0];
sx q[0];
rz(1.0157981) q[0];
rz(-0.47758684) q[2];
sx q[2];
rz(-1.4147593) q[2];
sx q[2];
rz(2.0618771) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.7970711) q[1];
sx q[1];
rz(-0.34725571) q[1];
sx q[1];
rz(-3.0879435) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.82181886) q[3];
sx q[3];
rz(-2.3082409) q[3];
sx q[3];
rz(2.5488473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3064208) q[2];
sx q[2];
rz(-2.1495337) q[2];
sx q[2];
rz(-1.4527808) q[2];
rz(0.49550223) q[3];
sx q[3];
rz(-0.82363868) q[3];
sx q[3];
rz(2.5775094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.751916) q[0];
sx q[0];
rz(-0.037411995) q[0];
sx q[0];
rz(2.9801242) q[0];
rz(-3.1096733) q[1];
sx q[1];
rz(-0.64161623) q[1];
sx q[1];
rz(-1.2772824) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.204435) q[0];
sx q[0];
rz(-1.8291408) q[0];
sx q[0];
rz(0.72465557) q[0];
rz(-pi) q[1];
rz(1.1487537) q[2];
sx q[2];
rz(-1.7817871) q[2];
sx q[2];
rz(0.27560292) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.016170597) q[1];
sx q[1];
rz(-1.6787663) q[1];
sx q[1];
rz(-2.8362464) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6971436) q[3];
sx q[3];
rz(-2.845361) q[3];
sx q[3];
rz(1.1100696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6892467) q[2];
sx q[2];
rz(-0.9114868) q[2];
sx q[2];
rz(-2.7455043) q[2];
rz(-0.39997697) q[3];
sx q[3];
rz(-0.59498274) q[3];
sx q[3];
rz(2.2649435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24092995) q[0];
sx q[0];
rz(-3.1234968) q[0];
sx q[0];
rz(-2.990429) q[0];
rz(0.97686544) q[1];
sx q[1];
rz(-0.38115373) q[1];
sx q[1];
rz(-0.74473286) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3522986) q[0];
sx q[0];
rz(-1.2937102) q[0];
sx q[0];
rz(-2.9267071) q[0];
rz(-1.2538386) q[2];
sx q[2];
rz(-0.43192568) q[2];
sx q[2];
rz(-3.0639086) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3724461) q[1];
sx q[1];
rz(-2.219104) q[1];
sx q[1];
rz(-1.5689956) q[1];
x q[2];
rz(2.7895097) q[3];
sx q[3];
rz(-1.7976645) q[3];
sx q[3];
rz(2.9405103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.62290827) q[2];
sx q[2];
rz(-1.9025977) q[2];
sx q[2];
rz(2.1935479) q[2];
rz(-1.0575804) q[3];
sx q[3];
rz(-1.6653929) q[3];
sx q[3];
rz(2.0675366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.093216151) q[0];
sx q[0];
rz(-0.26234782) q[0];
sx q[0];
rz(-0.23396215) q[0];
rz(-1.9562862) q[1];
sx q[1];
rz(-2.2133841) q[1];
sx q[1];
rz(-1.8918461) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5025185) q[0];
sx q[0];
rz(-0.6393798) q[0];
sx q[0];
rz(2.6431497) q[0];
rz(-pi) q[1];
x q[1];
rz(0.62057497) q[2];
sx q[2];
rz(-1.7399551) q[2];
sx q[2];
rz(-2.754654) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0873336) q[1];
sx q[1];
rz(-2.1076084) q[1];
sx q[1];
rz(-0.76294239) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.21548157) q[3];
sx q[3];
rz(-2.4487447) q[3];
sx q[3];
rz(2.1699435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.027111) q[2];
sx q[2];
rz(-1.9805084) q[2];
sx q[2];
rz(-1.4410045) q[2];
rz(0.13127413) q[3];
sx q[3];
rz(-2.6797397) q[3];
sx q[3];
rz(-2.89768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.048653614) q[0];
sx q[0];
rz(-0.96554148) q[0];
sx q[0];
rz(1.1337793) q[0];
rz(-2.1009905) q[1];
sx q[1];
rz(-0.84747172) q[1];
sx q[1];
rz(-0.52451959) q[1];
rz(-1.16083) q[2];
sx q[2];
rz(-0.66902918) q[2];
sx q[2];
rz(1.4794028) q[2];
rz(2.9374801) q[3];
sx q[3];
rz(-0.64328803) q[3];
sx q[3];
rz(-2.358123) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
