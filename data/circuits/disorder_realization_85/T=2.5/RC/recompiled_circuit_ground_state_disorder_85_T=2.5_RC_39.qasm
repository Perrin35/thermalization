OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.50897941) q[0];
sx q[0];
rz(-0.84492961) q[0];
sx q[0];
rz(0.71339575) q[0];
rz(-0.20819918) q[1];
sx q[1];
rz(-2.9543076) q[1];
sx q[1];
rz(2.0050144) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5644585) q[0];
sx q[0];
rz(-2.4043879) q[0];
sx q[0];
rz(-1.9335102) q[0];
rz(2.2520653) q[2];
sx q[2];
rz(-1.5233381) q[2];
sx q[2];
rz(1.0575546) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9507192) q[1];
sx q[1];
rz(-1.356186) q[1];
sx q[1];
rz(1.4282272) q[1];
x q[2];
rz(-1.7039812) q[3];
sx q[3];
rz(-1.2696506) q[3];
sx q[3];
rz(-0.74046053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1897757) q[2];
sx q[2];
rz(-1.8686029) q[2];
sx q[2];
rz(-0.58958685) q[2];
rz(-0.014852614) q[3];
sx q[3];
rz(-1.2735406) q[3];
sx q[3];
rz(-2.5428037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0066765229) q[0];
sx q[0];
rz(-2.0966457) q[0];
sx q[0];
rz(2.6892804) q[0];
rz(-2.2882838) q[1];
sx q[1];
rz(-2.6494458) q[1];
sx q[1];
rz(-2.7795627) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89328448) q[0];
sx q[0];
rz(-2.4442857) q[0];
sx q[0];
rz(2.6261283) q[0];
x q[1];
rz(0.49136038) q[2];
sx q[2];
rz(-1.7202913) q[2];
sx q[2];
rz(-2.4765365) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5278261) q[1];
sx q[1];
rz(-0.98694618) q[1];
sx q[1];
rz(0.85090249) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.70254247) q[3];
sx q[3];
rz(-1.6857393) q[3];
sx q[3];
rz(-3.0963932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8028458) q[2];
sx q[2];
rz(-0.57463988) q[2];
sx q[2];
rz(-0.90633264) q[2];
rz(-1.4551506) q[3];
sx q[3];
rz(-2.1833799) q[3];
sx q[3];
rz(1.5549392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0543268) q[0];
sx q[0];
rz(-1.9540906) q[0];
sx q[0];
rz(-2.1677256) q[0];
rz(-2.5663238) q[1];
sx q[1];
rz(-1.560805) q[1];
sx q[1];
rz(1.5781933) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3590976) q[0];
sx q[0];
rz(-1.5486915) q[0];
sx q[0];
rz(1.5074411) q[0];
rz(-pi) q[1];
rz(2.2996284) q[2];
sx q[2];
rz(-1.2281111) q[2];
sx q[2];
rz(-0.9093284) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8969054) q[1];
sx q[1];
rz(-0.88996038) q[1];
sx q[1];
rz(-0.61086706) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9138278) q[3];
sx q[3];
rz(-2.357759) q[3];
sx q[3];
rz(-1.7826102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.32597184) q[2];
sx q[2];
rz(-1.9591363) q[2];
sx q[2];
rz(-1.6131442) q[2];
rz(-0.01550393) q[3];
sx q[3];
rz(-1.9663234) q[3];
sx q[3];
rz(1.4874682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4267321) q[0];
sx q[0];
rz(-2.4437014) q[0];
sx q[0];
rz(-0.062407169) q[0];
rz(-1.1563673) q[1];
sx q[1];
rz(-0.9461177) q[1];
sx q[1];
rz(0.83414331) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7251563) q[0];
sx q[0];
rz(-1.8653989) q[0];
sx q[0];
rz(-1.0480773) q[0];
rz(1.4445199) q[2];
sx q[2];
rz(-1.6202721) q[2];
sx q[2];
rz(2.6425608) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4247113) q[1];
sx q[1];
rz(-2.5411027) q[1];
sx q[1];
rz(1.475715) q[1];
rz(-1.2150202) q[3];
sx q[3];
rz(-0.54635433) q[3];
sx q[3];
rz(2.234573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1674898) q[2];
sx q[2];
rz(-0.17593273) q[2];
sx q[2];
rz(-1.2220194) q[2];
rz(-0.40056285) q[3];
sx q[3];
rz(-2.1192854) q[3];
sx q[3];
rz(0.21489531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.062531384) q[0];
sx q[0];
rz(-2.489594) q[0];
sx q[0];
rz(-2.8849211) q[0];
rz(-1.7968862) q[1];
sx q[1];
rz(-1.2213629) q[1];
sx q[1];
rz(1.409097) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.870011) q[0];
sx q[0];
rz(-0.74423941) q[0];
sx q[0];
rz(1.1209473) q[0];
rz(2.1878384) q[2];
sx q[2];
rz(-0.34025345) q[2];
sx q[2];
rz(-2.798852) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0875279) q[1];
sx q[1];
rz(-1.9913186) q[1];
sx q[1];
rz(-1.5271212) q[1];
rz(-1.2692637) q[3];
sx q[3];
rz(-0.9253987) q[3];
sx q[3];
rz(-2.0567578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.04074374) q[2];
sx q[2];
rz(-0.70793968) q[2];
sx q[2];
rz(2.9948998) q[2];
rz(1.7723068) q[3];
sx q[3];
rz(-2.7761603) q[3];
sx q[3];
rz(-0.2379612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.20928243) q[0];
sx q[0];
rz(-0.76618659) q[0];
sx q[0];
rz(0.76195088) q[0];
rz(-0.85488287) q[1];
sx q[1];
rz(-1.8652752) q[1];
sx q[1];
rz(-2.4529822) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.541674) q[0];
sx q[0];
rz(-1.1038938) q[0];
sx q[0];
rz(3.0476157) q[0];
rz(-2.9607331) q[2];
sx q[2];
rz(-3.0263889) q[2];
sx q[2];
rz(1.2173139) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2999867) q[1];
sx q[1];
rz(-3.0191023) q[1];
sx q[1];
rz(2.5754506) q[1];
rz(0.33276673) q[3];
sx q[3];
rz(-1.1586572) q[3];
sx q[3];
rz(-1.0650959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.40053549) q[2];
sx q[2];
rz(-2.1964938) q[2];
sx q[2];
rz(0.98089027) q[2];
rz(0.26997057) q[3];
sx q[3];
rz(-1.910784) q[3];
sx q[3];
rz(0.009416906) q[3];
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
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16801676) q[0];
sx q[0];
rz(-1.4485757) q[0];
sx q[0];
rz(0.5434522) q[0];
rz(1.5914397) q[1];
sx q[1];
rz(-0.88528577) q[1];
sx q[1];
rz(-0.98522225) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7783981) q[0];
sx q[0];
rz(-1.8002602) q[0];
sx q[0];
rz(2.6465332) q[0];
rz(-2.1210873) q[2];
sx q[2];
rz(-0.63945635) q[2];
sx q[2];
rz(2.8957862) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0646202) q[1];
sx q[1];
rz(-2.1482921) q[1];
sx q[1];
rz(-2.6813238) q[1];
rz(-pi) q[2];
rz(0.51237088) q[3];
sx q[3];
rz(-1.8647883) q[3];
sx q[3];
rz(-1.1666573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8103509) q[2];
sx q[2];
rz(-2.6369429) q[2];
sx q[2];
rz(-1.3387559) q[2];
rz(1.6797558) q[3];
sx q[3];
rz(-1.5509408) q[3];
sx q[3];
rz(-0.67080474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93002334) q[0];
sx q[0];
rz(-1.1299364) q[0];
sx q[0];
rz(1.2343963) q[0];
rz(-1.3537539) q[1];
sx q[1];
rz(-2.6852971) q[1];
sx q[1];
rz(-3.0455468) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98063606) q[0];
sx q[0];
rz(-0.90117878) q[0];
sx q[0];
rz(-2.5494844) q[0];
x q[1];
rz(-3.0158218) q[2];
sx q[2];
rz(-0.54361225) q[2];
sx q[2];
rz(2.8961492) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7754566) q[1];
sx q[1];
rz(-1.6153496) q[1];
sx q[1];
rz(0.38031995) q[1];
rz(-pi) q[2];
rz(0.84026321) q[3];
sx q[3];
rz(-1.6018768) q[3];
sx q[3];
rz(-0.66339359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.3357521) q[2];
sx q[2];
rz(-0.38613191) q[2];
sx q[2];
rz(-1.1248355) q[2];
rz(-0.8791033) q[3];
sx q[3];
rz(-1.2259038) q[3];
sx q[3];
rz(0.44900289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45824555) q[0];
sx q[0];
rz(-1.6793716) q[0];
sx q[0];
rz(0.71518389) q[0];
rz(0.31582754) q[1];
sx q[1];
rz(-0.68604699) q[1];
sx q[1];
rz(-0.16206965) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.049555) q[0];
sx q[0];
rz(-0.83956912) q[0];
sx q[0];
rz(-0.58734307) q[0];
rz(-pi) q[1];
rz(-2.8388073) q[2];
sx q[2];
rz(-2.5368779) q[2];
sx q[2];
rz(2.5446937) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4906076) q[1];
sx q[1];
rz(-0.7443634) q[1];
sx q[1];
rz(2.8815844) q[1];
rz(-pi) q[2];
rz(-1.9088184) q[3];
sx q[3];
rz(-1.2673339) q[3];
sx q[3];
rz(0.4952411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7519303) q[2];
sx q[2];
rz(-1.1952362) q[2];
sx q[2];
rz(1.2644838) q[2];
rz(1.3324995) q[3];
sx q[3];
rz(-1.8615362) q[3];
sx q[3];
rz(-0.31681028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7952809) q[0];
sx q[0];
rz(-0.77487159) q[0];
sx q[0];
rz(1.0493904) q[0];
rz(-0.28114444) q[1];
sx q[1];
rz(-2.099359) q[1];
sx q[1];
rz(-1.8195456) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3157177) q[0];
sx q[0];
rz(-1.0260887) q[0];
sx q[0];
rz(1.8648575) q[0];
rz(-2.5192267) q[2];
sx q[2];
rz(-0.73957788) q[2];
sx q[2];
rz(-0.057829521) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2019649) q[1];
sx q[1];
rz(-1.2117821) q[1];
sx q[1];
rz(-0.16965349) q[1];
rz(2.9024944) q[3];
sx q[3];
rz(-0.39632583) q[3];
sx q[3];
rz(-2.9843085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.37810024) q[2];
sx q[2];
rz(-1.3214448) q[2];
sx q[2];
rz(-3.0901001) q[2];
rz(1.817305) q[3];
sx q[3];
rz(-1.6591502) q[3];
sx q[3];
rz(-1.2073368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0180203) q[0];
sx q[0];
rz(-1.3326895) q[0];
sx q[0];
rz(1.4154758) q[0];
rz(-0.82072683) q[1];
sx q[1];
rz(-1.6976994) q[1];
sx q[1];
rz(1.2669947) q[1];
rz(-2.5222798) q[2];
sx q[2];
rz(-1.3765341) q[2];
sx q[2];
rz(1.7159009) q[2];
rz(0.5140082) q[3];
sx q[3];
rz(-1.2759943) q[3];
sx q[3];
rz(-1.3392824) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
