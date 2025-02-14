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
rz(1.9198298) q[0];
sx q[0];
rz(-1.3718995) q[0];
sx q[0];
rz(-2.0576117) q[0];
rz(2.3594175) q[1];
sx q[1];
rz(-0.44936925) q[1];
sx q[1];
rz(2.3545797) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0660853) q[0];
sx q[0];
rz(-1.721816) q[0];
sx q[0];
rz(1.9635033) q[0];
rz(-pi) q[1];
rz(-1.6264362) q[2];
sx q[2];
rz(-0.67650992) q[2];
sx q[2];
rz(1.8970053) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2725153) q[1];
sx q[1];
rz(-2.1109588) q[1];
sx q[1];
rz(-0.16611986) q[1];
x q[2];
rz(2.8008599) q[3];
sx q[3];
rz(-2.2267146) q[3];
sx q[3];
rz(-0.64217303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3419753) q[2];
sx q[2];
rz(-1.3222597) q[2];
sx q[2];
rz(2.4413617) q[2];
rz(1.7796984) q[3];
sx q[3];
rz(-1.9340065) q[3];
sx q[3];
rz(3.0119058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48548651) q[0];
sx q[0];
rz(-0.34657297) q[0];
sx q[0];
rz(-1.9488126) q[0];
rz(2.9876409) q[1];
sx q[1];
rz(-0.62869453) q[1];
sx q[1];
rz(-1.2492294) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9500268) q[0];
sx q[0];
rz(-0.18432549) q[0];
sx q[0];
rz(-0.34785716) q[0];
rz(-pi) q[1];
x q[1];
rz(0.05441101) q[2];
sx q[2];
rz(-2.2073032) q[2];
sx q[2];
rz(2.7639703) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2613879) q[1];
sx q[1];
rz(-1.4900165) q[1];
sx q[1];
rz(-2.8360663) q[1];
x q[2];
rz(2.0286247) q[3];
sx q[3];
rz(-2.4220123) q[3];
sx q[3];
rz(-0.45872575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2092756) q[2];
sx q[2];
rz(-1.2195769) q[2];
sx q[2];
rz(0.90927124) q[2];
rz(-0.13898177) q[3];
sx q[3];
rz(-2.8863972) q[3];
sx q[3];
rz(-2.2709258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22055498) q[0];
sx q[0];
rz(-1.166936) q[0];
sx q[0];
rz(-1.9535109) q[0];
rz(-2.2782169) q[1];
sx q[1];
rz(-1.8977576) q[1];
sx q[1];
rz(-0.9695425) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70308951) q[0];
sx q[0];
rz(-2.9155779) q[0];
sx q[0];
rz(-0.78031197) q[0];
rz(-pi) q[1];
rz(1.9447504) q[2];
sx q[2];
rz(-1.8965169) q[2];
sx q[2];
rz(2.0545127) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4789106) q[1];
sx q[1];
rz(-1.3363774) q[1];
sx q[1];
rz(1.6699416) q[1];
x q[2];
rz(2.5695877) q[3];
sx q[3];
rz(-2.1915132) q[3];
sx q[3];
rz(-0.86670029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.36677507) q[2];
sx q[2];
rz(-0.25664169) q[2];
sx q[2];
rz(2.4859264) q[2];
rz(-1.5244779) q[3];
sx q[3];
rz(-1.3527801) q[3];
sx q[3];
rz(-0.85339439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.26576385) q[0];
sx q[0];
rz(-2.0142856) q[0];
sx q[0];
rz(1.6831336) q[0];
rz(1.4153076) q[1];
sx q[1];
rz(-1.4348607) q[1];
sx q[1];
rz(1.2687792) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.388871) q[0];
sx q[0];
rz(-1.9170887) q[0];
sx q[0];
rz(-0.66360051) q[0];
rz(-pi) q[1];
x q[1];
rz(0.54019467) q[2];
sx q[2];
rz(-1.0120076) q[2];
sx q[2];
rz(-1.9524198) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.1139196) q[1];
sx q[1];
rz(-2.259127) q[1];
sx q[1];
rz(2.9034241) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7990718) q[3];
sx q[3];
rz(-1.6581495) q[3];
sx q[3];
rz(-2.563208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.65712523) q[2];
sx q[2];
rz(-1.7871658) q[2];
sx q[2];
rz(1.2628868) q[2];
rz(-2.9386988) q[3];
sx q[3];
rz(-1.032136) q[3];
sx q[3];
rz(-1.8111022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48388594) q[0];
sx q[0];
rz(-1.0840451) q[0];
sx q[0];
rz(0.44802353) q[0];
rz(1.9810642) q[1];
sx q[1];
rz(-1.0803761) q[1];
sx q[1];
rz(-2.0652658) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3884013) q[0];
sx q[0];
rz(-1.7832169) q[0];
sx q[0];
rz(2.302049) q[0];
rz(0.33319039) q[2];
sx q[2];
rz(-1.6181862) q[2];
sx q[2];
rz(-1.2409925) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.68491918) q[1];
sx q[1];
rz(-0.071487814) q[1];
sx q[1];
rz(-1.0409357) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5994051) q[3];
sx q[3];
rz(-2.5243658) q[3];
sx q[3];
rz(2.9742179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6059604) q[2];
sx q[2];
rz(-1.2781906) q[2];
sx q[2];
rz(2.7995301) q[2];
rz(-2.0096807) q[3];
sx q[3];
rz(-1.071238) q[3];
sx q[3];
rz(-0.38558495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.015942052) q[0];
sx q[0];
rz(-2.5378939) q[0];
sx q[0];
rz(-1.1018671) q[0];
rz(2.8562538) q[1];
sx q[1];
rz(-0.97831786) q[1];
sx q[1];
rz(-2.2314821) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9302437) q[0];
sx q[0];
rz(-1.6841297) q[0];
sx q[0];
rz(-1.428753) q[0];
x q[1];
rz(1.455972) q[2];
sx q[2];
rz(-1.6636724) q[2];
sx q[2];
rz(0.50485134) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8273298) q[1];
sx q[1];
rz(-2.0138965) q[1];
sx q[1];
rz(3.0783975) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7496787) q[3];
sx q[3];
rz(-1.4174882) q[3];
sx q[3];
rz(1.2962504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.82112271) q[2];
sx q[2];
rz(-1.4218825) q[2];
sx q[2];
rz(1.0538496) q[2];
rz(0.36367917) q[3];
sx q[3];
rz(-1.4854919) q[3];
sx q[3];
rz(0.53500879) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66154552) q[0];
sx q[0];
rz(-1.0129901) q[0];
sx q[0];
rz(0.15740982) q[0];
rz(-0.62955457) q[1];
sx q[1];
rz(-1.5512356) q[1];
sx q[1];
rz(-0.079364337) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4422121) q[0];
sx q[0];
rz(-1.4794425) q[0];
sx q[0];
rz(-0.070789074) q[0];
rz(-pi) q[1];
rz(-2.7474097) q[2];
sx q[2];
rz(-1.7346343) q[2];
sx q[2];
rz(-0.25111408) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.094202894) q[1];
sx q[1];
rz(-1.6999287) q[1];
sx q[1];
rz(-2.6455621) q[1];
rz(1.5299876) q[3];
sx q[3];
rz(-2.5548078) q[3];
sx q[3];
rz(-0.87848488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.51084149) q[2];
sx q[2];
rz(-1.3151104) q[2];
sx q[2];
rz(-0.17244478) q[2];
rz(2.4991961) q[3];
sx q[3];
rz(-0.96697092) q[3];
sx q[3];
rz(2.7914458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.073535) q[0];
sx q[0];
rz(-1.8386766) q[0];
sx q[0];
rz(0.56655836) q[0];
rz(2.9134275) q[1];
sx q[1];
rz(-1.4527495) q[1];
sx q[1];
rz(-1.8295005) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8926548) q[0];
sx q[0];
rz(-2.3088107) q[0];
sx q[0];
rz(-2.8036462) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.84374441) q[2];
sx q[2];
rz(-1.4829501) q[2];
sx q[2];
rz(1.530046) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6375288) q[1];
sx q[1];
rz(-2.3603068) q[1];
sx q[1];
rz(-2.3243796) q[1];
x q[2];
rz(3.1080179) q[3];
sx q[3];
rz(-1.9596425) q[3];
sx q[3];
rz(-1.7186708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.18099004) q[2];
sx q[2];
rz(-1.484551) q[2];
sx q[2];
rz(-2.9385938) q[2];
rz(2.3024043) q[3];
sx q[3];
rz(-2.8036717) q[3];
sx q[3];
rz(2.2904229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7508271) q[0];
sx q[0];
rz(-0.44742328) q[0];
sx q[0];
rz(0.14697337) q[0];
rz(-0.43288639) q[1];
sx q[1];
rz(-1.9501481) q[1];
sx q[1];
rz(-1.254522) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1330652) q[0];
sx q[0];
rz(-2.0957895) q[0];
sx q[0];
rz(1.1675444) q[0];
rz(-pi) q[1];
rz(-1.8190838) q[2];
sx q[2];
rz(-1.5380368) q[2];
sx q[2];
rz(-2.0078307) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.93941037) q[1];
sx q[1];
rz(-1.2263884) q[1];
sx q[1];
rz(2.9581974) q[1];
rz(-pi) q[2];
rz(1.5490554) q[3];
sx q[3];
rz(-1.4026902) q[3];
sx q[3];
rz(-2.4631207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7811084) q[2];
sx q[2];
rz(-1.0177178) q[2];
sx q[2];
rz(0.54070365) q[2];
rz(0.30132076) q[3];
sx q[3];
rz(-0.40073985) q[3];
sx q[3];
rz(1.8436684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28894579) q[0];
sx q[0];
rz(-1.1971373) q[0];
sx q[0];
rz(0.72809118) q[0];
rz(-2.6241265) q[1];
sx q[1];
rz(-1.4897852) q[1];
sx q[1];
rz(2.1450086) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40328011) q[0];
sx q[0];
rz(-1.9150253) q[0];
sx q[0];
rz(-3.1136373) q[0];
x q[1];
rz(-1.4696737) q[2];
sx q[2];
rz(-0.99144672) q[2];
sx q[2];
rz(0.92321009) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0658465) q[1];
sx q[1];
rz(-1.5567413) q[1];
sx q[1];
rz(-1.3142287) q[1];
rz(1.4703937) q[3];
sx q[3];
rz(-1.3884423) q[3];
sx q[3];
rz(-0.34913464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1835798) q[2];
sx q[2];
rz(-1.5637584) q[2];
sx q[2];
rz(0.070240423) q[2];
rz(0.0066512935) q[3];
sx q[3];
rz(-0.17156048) q[3];
sx q[3];
rz(0.10449617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2086647) q[0];
sx q[0];
rz(-1.7345971) q[0];
sx q[0];
rz(-1.4165184) q[0];
rz(2.3729462) q[1];
sx q[1];
rz(-1.9192764) q[1];
sx q[1];
rz(2.8142014) q[1];
rz(-1.6835536) q[2];
sx q[2];
rz(-1.590645) q[2];
sx q[2];
rz(-2.0525862) q[2];
rz(-0.27003084) q[3];
sx q[3];
rz(-1.217801) q[3];
sx q[3];
rz(1.4767811) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
