OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7497082) q[0];
sx q[0];
rz(-2.9449129) q[0];
sx q[0];
rz(-1.1893907) q[0];
rz(-2.9130798) q[1];
sx q[1];
rz(-2.300188) q[1];
sx q[1];
rz(2.7639311) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1188287) q[0];
sx q[0];
rz(-1.6435992) q[0];
sx q[0];
rz(1.692549) q[0];
rz(-1.471465) q[2];
sx q[2];
rz(-2.8566395) q[2];
sx q[2];
rz(0.003665912) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.764896) q[1];
sx q[1];
rz(-2.2590859) q[1];
sx q[1];
rz(-2.1535758) q[1];
rz(2.8336485) q[3];
sx q[3];
rz(-0.43711284) q[3];
sx q[3];
rz(-1.8969769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1315786) q[2];
sx q[2];
rz(-0.49396124) q[2];
sx q[2];
rz(-0.67260355) q[2];
rz(0.16942313) q[3];
sx q[3];
rz(-0.38893458) q[3];
sx q[3];
rz(1.3385564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23068962) q[0];
sx q[0];
rz(-2.2407273) q[0];
sx q[0];
rz(2.9130274) q[0];
rz(2.9810492) q[1];
sx q[1];
rz(-1.4385782) q[1];
sx q[1];
rz(0.28796089) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1624958) q[0];
sx q[0];
rz(-1.393035) q[0];
sx q[0];
rz(1.7642679) q[0];
rz(-pi) q[1];
rz(2.7520913) q[2];
sx q[2];
rz(-2.1283538) q[2];
sx q[2];
rz(0.26706375) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8781232) q[1];
sx q[1];
rz(-0.97349226) q[1];
sx q[1];
rz(1.0695446) q[1];
rz(-1.6091299) q[3];
sx q[3];
rz(-1.2206856) q[3];
sx q[3];
rz(2.4083174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5471389) q[2];
sx q[2];
rz(-1.889166) q[2];
sx q[2];
rz(-1.4734369) q[2];
rz(-0.93747059) q[3];
sx q[3];
rz(-2.6963186) q[3];
sx q[3];
rz(-0.37500769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.32524747) q[0];
sx q[0];
rz(-0.49961093) q[0];
sx q[0];
rz(0.26741272) q[0];
rz(1.7193517) q[1];
sx q[1];
rz(-1.1317252) q[1];
sx q[1];
rz(-0.95169383) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50599498) q[0];
sx q[0];
rz(-1.4800758) q[0];
sx q[0];
rz(1.7631625) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9159413) q[2];
sx q[2];
rz(-1.0661085) q[2];
sx q[2];
rz(0.80863189) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.051085) q[1];
sx q[1];
rz(-1.8790434) q[1];
sx q[1];
rz(-2.3329263) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1943201) q[3];
sx q[3];
rz(-1.1215278) q[3];
sx q[3];
rz(2.505213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0114228) q[2];
sx q[2];
rz(-1.495785) q[2];
sx q[2];
rz(0.34417957) q[2];
rz(2.2551645) q[3];
sx q[3];
rz(-2.8635946) q[3];
sx q[3];
rz(-2.5755431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13720559) q[0];
sx q[0];
rz(-1.6442278) q[0];
sx q[0];
rz(1.8970867) q[0];
rz(0.1291153) q[1];
sx q[1];
rz(-1.3543509) q[1];
sx q[1];
rz(2.7688162) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6308206) q[0];
sx q[0];
rz(-1.456506) q[0];
sx q[0];
rz(0.75829102) q[0];
x q[1];
rz(-0.77809019) q[2];
sx q[2];
rz(-1.8847244) q[2];
sx q[2];
rz(0.84184605) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7074234) q[1];
sx q[1];
rz(-1.3278828) q[1];
sx q[1];
rz(1.0711819) q[1];
x q[2];
rz(-0.32897207) q[3];
sx q[3];
rz(-1.1165459) q[3];
sx q[3];
rz(0.3262375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.63561511) q[2];
sx q[2];
rz(-1.4900692) q[2];
sx q[2];
rz(0.90488952) q[2];
rz(-2.7010226) q[3];
sx q[3];
rz(-2.6456656) q[3];
sx q[3];
rz(-2.1499965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6624517) q[0];
sx q[0];
rz(-2.9859556) q[0];
sx q[0];
rz(1.8537846) q[0];
rz(1.4783391) q[1];
sx q[1];
rz(-1.1454502) q[1];
sx q[1];
rz(-0.53422654) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68191093) q[0];
sx q[0];
rz(-0.43897438) q[0];
sx q[0];
rz(0.79553332) q[0];
rz(-pi) q[1];
rz(-2.1557501) q[2];
sx q[2];
rz(-2.5702229) q[2];
sx q[2];
rz(2.5845598) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5200107) q[1];
sx q[1];
rz(-0.098712155) q[1];
sx q[1];
rz(2.2732544) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4228135) q[3];
sx q[3];
rz(-1.4712417) q[3];
sx q[3];
rz(2.0841141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6405032) q[2];
sx q[2];
rz(-1.2155632) q[2];
sx q[2];
rz(0.14349288) q[2];
rz(1.3714553) q[3];
sx q[3];
rz(-1.8618795) q[3];
sx q[3];
rz(0.21970704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4390398) q[0];
sx q[0];
rz(-0.87600231) q[0];
sx q[0];
rz(-0.23705661) q[0];
rz(1.2409695) q[1];
sx q[1];
rz(-0.82890141) q[1];
sx q[1];
rz(-1.7664849) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28778827) q[0];
sx q[0];
rz(-1.8332991) q[0];
sx q[0];
rz(2.7909423) q[0];
x q[1];
rz(1.6429971) q[2];
sx q[2];
rz(-2.318396) q[2];
sx q[2];
rz(0.88897926) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4543537) q[1];
sx q[1];
rz(-0.41077405) q[1];
sx q[1];
rz(-0.06128581) q[1];
x q[2];
rz(1.1477908) q[3];
sx q[3];
rz(-2.0241963) q[3];
sx q[3];
rz(-0.51844937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.53283006) q[2];
sx q[2];
rz(-2.4217748) q[2];
sx q[2];
rz(0.14870816) q[2];
rz(-0.016629774) q[3];
sx q[3];
rz(-0.36706585) q[3];
sx q[3];
rz(-0.19255157) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49611133) q[0];
sx q[0];
rz(-2.8493024) q[0];
sx q[0];
rz(-2.2684229) q[0];
rz(-2.219615) q[1];
sx q[1];
rz(-2.0261804) q[1];
sx q[1];
rz(-3.057664) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.98691) q[0];
sx q[0];
rz(-2.1213518) q[0];
sx q[0];
rz(-0.26225342) q[0];
rz(-pi) q[1];
rz(1.3533808) q[2];
sx q[2];
rz(-1.916421) q[2];
sx q[2];
rz(-2.2923922) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1740239) q[1];
sx q[1];
rz(-2.1597383) q[1];
sx q[1];
rz(0.80935652) q[1];
rz(-pi) q[2];
rz(0.84079068) q[3];
sx q[3];
rz(-1.519324) q[3];
sx q[3];
rz(-1.7712902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7579047) q[2];
sx q[2];
rz(-2.7010475) q[2];
sx q[2];
rz(2.596358) q[2];
rz(-2.7111354) q[3];
sx q[3];
rz(-1.1377708) q[3];
sx q[3];
rz(2.8366413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51171821) q[0];
sx q[0];
rz(-0.015462333) q[0];
sx q[0];
rz(-1.2114725) q[0];
rz(2.1684872) q[1];
sx q[1];
rz(-2.548023) q[1];
sx q[1];
rz(2.1957695) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5686533) q[0];
sx q[0];
rz(-1.6488254) q[0];
sx q[0];
rz(-0.22139876) q[0];
rz(-pi) q[1];
rz(2.4589355) q[2];
sx q[2];
rz(-0.97634146) q[2];
sx q[2];
rz(0.97460954) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.91499) q[1];
sx q[1];
rz(-1.9493002) q[1];
sx q[1];
rz(-1.539991) q[1];
rz(0.66611992) q[3];
sx q[3];
rz(-1.564581) q[3];
sx q[3];
rz(-0.70762779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2716081) q[2];
sx q[2];
rz(-1.4166069) q[2];
sx q[2];
rz(-2.8611709) q[2];
rz(0.60493207) q[3];
sx q[3];
rz(-1.0792024) q[3];
sx q[3];
rz(0.98208565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1542926) q[0];
sx q[0];
rz(-0.91006088) q[0];
sx q[0];
rz(-2.5262685) q[0];
rz(-0.92957169) q[1];
sx q[1];
rz(-0.85837448) q[1];
sx q[1];
rz(-2.5659134) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.070411365) q[0];
sx q[0];
rz(-2.1343263) q[0];
sx q[0];
rz(2.1348743) q[0];
x q[1];
rz(0.11328463) q[2];
sx q[2];
rz(-1.3759398) q[2];
sx q[2];
rz(1.4154797) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.5902192) q[1];
sx q[1];
rz(-3.0312523) q[1];
sx q[1];
rz(-1.3361841) q[1];
rz(1.0300893) q[3];
sx q[3];
rz(-1.732015) q[3];
sx q[3];
rz(0.89564039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3118887) q[2];
sx q[2];
rz(-0.76278937) q[2];
sx q[2];
rz(-3.1196307) q[2];
rz(2.9711376) q[3];
sx q[3];
rz(-2.1178092) q[3];
sx q[3];
rz(-0.26486614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4158674) q[0];
sx q[0];
rz(-0.9265582) q[0];
sx q[0];
rz(2.5073994) q[0];
rz(3.1126853) q[1];
sx q[1];
rz(-0.78866619) q[1];
sx q[1];
rz(2.172487) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.076981912) q[0];
sx q[0];
rz(-1.4990028) q[0];
sx q[0];
rz(-0.063477091) q[0];
x q[1];
rz(-0.74659851) q[2];
sx q[2];
rz(-1.4274297) q[2];
sx q[2];
rz(1.5310841) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9173968) q[1];
sx q[1];
rz(-2.2317413) q[1];
sx q[1];
rz(-1.7978653) q[1];
rz(-pi) q[2];
rz(2.7006847) q[3];
sx q[3];
rz(-1.3135859) q[3];
sx q[3];
rz(-2.5901026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6897631) q[2];
sx q[2];
rz(-1.657106) q[2];
sx q[2];
rz(0.36007145) q[2];
rz(1.6137971) q[3];
sx q[3];
rz(-2.4751622) q[3];
sx q[3];
rz(2.578919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9344899) q[0];
sx q[0];
rz(-1.5710545) q[0];
sx q[0];
rz(1.5221773) q[0];
rz(3.0974401) q[1];
sx q[1];
rz(-1.4587198) q[1];
sx q[1];
rz(-1.1062467) q[1];
rz(0.85514851) q[2];
sx q[2];
rz(-2.6519041) q[2];
sx q[2];
rz(2.3408163) q[2];
rz(0.58892693) q[3];
sx q[3];
rz(-0.52994655) q[3];
sx q[3];
rz(0.51701057) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
