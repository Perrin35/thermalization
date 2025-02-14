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
rz(0.084963381) q[0];
sx q[0];
rz(-2.8391916) q[0];
sx q[0];
rz(3.1095355) q[0];
rz(3.1337466) q[1];
sx q[1];
rz(3.6454522) q[1];
sx q[1];
rz(7.03581) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1330452) q[0];
sx q[0];
rz(-0.38131443) q[0];
sx q[0];
rz(-0.56579907) q[0];
x q[1];
rz(-0.36280323) q[2];
sx q[2];
rz(-1.9224836) q[2];
sx q[2];
rz(-1.2376518) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.17129414) q[1];
sx q[1];
rz(-0.42810218) q[1];
sx q[1];
rz(0.64117214) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9273622) q[3];
sx q[3];
rz(-1.9211244) q[3];
sx q[3];
rz(0.080834724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.96693119) q[2];
sx q[2];
rz(-1.175468) q[2];
sx q[2];
rz(2.3215129) q[2];
rz(-0.44102937) q[3];
sx q[3];
rz(-2.0377772) q[3];
sx q[3];
rz(0.57344121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.012861982) q[0];
sx q[0];
rz(-0.91330376) q[0];
sx q[0];
rz(-2.7313857) q[0];
rz(-2.7606616) q[1];
sx q[1];
rz(-1.5313287) q[1];
sx q[1];
rz(1.7832696) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5986431) q[0];
sx q[0];
rz(-0.99587593) q[0];
sx q[0];
rz(-2.6680114) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.50308268) q[2];
sx q[2];
rz(-1.0840992) q[2];
sx q[2];
rz(-0.27931914) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1413026) q[1];
sx q[1];
rz(-1.4030255) q[1];
sx q[1];
rz(-3.0649158) q[1];
x q[2];
rz(1.3639327) q[3];
sx q[3];
rz(-2.4914722) q[3];
sx q[3];
rz(0.93322414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7831948) q[2];
sx q[2];
rz(-2.7326549) q[2];
sx q[2];
rz(0.50296339) q[2];
rz(1.2991692) q[3];
sx q[3];
rz(-2.3830569) q[3];
sx q[3];
rz(-0.23183307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7314887) q[0];
sx q[0];
rz(-2.6298611) q[0];
sx q[0];
rz(2.1657535) q[0];
rz(-0.93636912) q[1];
sx q[1];
rz(-1.5872637) q[1];
sx q[1];
rz(0.13793129) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5463211) q[0];
sx q[0];
rz(-1.0965276) q[0];
sx q[0];
rz(-2.7300937) q[0];
rz(2.7439026) q[2];
sx q[2];
rz(-1.4639336) q[2];
sx q[2];
rz(1.7025089) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.960818) q[1];
sx q[1];
rz(-1.7021693) q[1];
sx q[1];
rz(-2.2714494) q[1];
rz(-pi) q[2];
x q[2];
rz(0.46584399) q[3];
sx q[3];
rz(-1.3261822) q[3];
sx q[3];
rz(2.1810093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5522573) q[2];
sx q[2];
rz(-1.6103123) q[2];
sx q[2];
rz(-1.662558) q[2];
rz(-0.79834437) q[3];
sx q[3];
rz(-0.84404293) q[3];
sx q[3];
rz(0.020309694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.8524356) q[0];
sx q[0];
rz(-3.030179) q[0];
sx q[0];
rz(2.074746) q[0];
rz(-1.5486708) q[1];
sx q[1];
rz(-1.2499864) q[1];
sx q[1];
rz(-0.41935316) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33808483) q[0];
sx q[0];
rz(-0.28059059) q[0];
sx q[0];
rz(-2.1795033) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.1544756) q[2];
sx q[2];
rz(-2.6082268) q[2];
sx q[2];
rz(1.5325002) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4511299) q[1];
sx q[1];
rz(-1.9277384) q[1];
sx q[1];
rz(-3.073077) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0556424) q[3];
sx q[3];
rz(-1.2618229) q[3];
sx q[3];
rz(-3.0687208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.35820094) q[2];
sx q[2];
rz(-0.50938598) q[2];
sx q[2];
rz(1.5436714) q[2];
rz(-1.915043) q[3];
sx q[3];
rz(-1.9251325) q[3];
sx q[3];
rz(-1.8515057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.2655547) q[0];
sx q[0];
rz(-2.4198678) q[0];
sx q[0];
rz(0.75620404) q[0];
rz(1.2528231) q[1];
sx q[1];
rz(-2.2105261) q[1];
sx q[1];
rz(-1.4160215) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21391103) q[0];
sx q[0];
rz(-2.1577303) q[0];
sx q[0];
rz(-1.3089433) q[0];
x q[1];
rz(-1.6966693) q[2];
sx q[2];
rz(-1.9121661) q[2];
sx q[2];
rz(-1.204513) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.49633138) q[1];
sx q[1];
rz(-2.1095536) q[1];
sx q[1];
rz(-1.5790198) q[1];
x q[2];
rz(2.7783423) q[3];
sx q[3];
rz(-1.7492069) q[3];
sx q[3];
rz(1.6476064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.34117928) q[2];
sx q[2];
rz(-0.74113733) q[2];
sx q[2];
rz(2.9608012) q[2];
rz(1.4888658) q[3];
sx q[3];
rz(-1.7427665) q[3];
sx q[3];
rz(1.6256049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7086696) q[0];
sx q[0];
rz(-2.0661418) q[0];
sx q[0];
rz(-2.3413626) q[0];
rz(2.8569787) q[1];
sx q[1];
rz(-2.0814643) q[1];
sx q[1];
rz(-0.12399331) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.894688) q[0];
sx q[0];
rz(-0.13291153) q[0];
sx q[0];
rz(-1.354643) q[0];
rz(-pi) q[1];
rz(1.8154549) q[2];
sx q[2];
rz(-2.8996433) q[2];
sx q[2];
rz(1.2173779) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.968144) q[1];
sx q[1];
rz(-0.3160797) q[1];
sx q[1];
rz(-2.8760002) q[1];
rz(-pi) q[2];
rz(-0.97449755) q[3];
sx q[3];
rz(-1.8829405) q[3];
sx q[3];
rz(-0.5398992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2064712) q[2];
sx q[2];
rz(-1.4149041) q[2];
sx q[2];
rz(0.07621152) q[2];
rz(1.291409) q[3];
sx q[3];
rz(-1.094787) q[3];
sx q[3];
rz(-0.61856234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9718219) q[0];
sx q[0];
rz(-1.562028) q[0];
sx q[0];
rz(-0.75702697) q[0];
rz(-0.79611671) q[1];
sx q[1];
rz(-1.4069822) q[1];
sx q[1];
rz(-0.98181358) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9011089) q[0];
sx q[0];
rz(-0.79190688) q[0];
sx q[0];
rz(-1.1674985) q[0];
rz(-pi) q[1];
rz(-3.0715838) q[2];
sx q[2];
rz(-1.4638136) q[2];
sx q[2];
rz(1.9221523) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.422324) q[1];
sx q[1];
rz(-2.0403721) q[1];
sx q[1];
rz(2.6828241) q[1];
rz(-pi) q[2];
x q[2];
rz(2.997274) q[3];
sx q[3];
rz(-1.8486406) q[3];
sx q[3];
rz(-1.2793737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.41165274) q[2];
sx q[2];
rz(-0.81092683) q[2];
sx q[2];
rz(-0.5438424) q[2];
rz(-0.020921556) q[3];
sx q[3];
rz(-1.436751) q[3];
sx q[3];
rz(-2.3232443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.2127317) q[0];
sx q[0];
rz(-0.91101557) q[0];
sx q[0];
rz(0.62335706) q[0];
rz(-2.8202672) q[1];
sx q[1];
rz(-0.61505452) q[1];
sx q[1];
rz(2.8772433) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67548114) q[0];
sx q[0];
rz(-2.3947885) q[0];
sx q[0];
rz(1.2561428) q[0];
rz(3.1000443) q[2];
sx q[2];
rz(-1.9061521) q[2];
sx q[2];
rz(1.2183508) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5874065) q[1];
sx q[1];
rz(-2.3465112) q[1];
sx q[1];
rz(1.9663747) q[1];
rz(-pi) q[2];
rz(-1.1846011) q[3];
sx q[3];
rz(-2.8004146) q[3];
sx q[3];
rz(2.6130849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.136772) q[2];
sx q[2];
rz(-2.4592168) q[2];
sx q[2];
rz(2.6668059) q[2];
rz(0.36561203) q[3];
sx q[3];
rz(-0.48706278) q[3];
sx q[3];
rz(0.18317187) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30597618) q[0];
sx q[0];
rz(-0.37537471) q[0];
sx q[0];
rz(-0.48582745) q[0];
rz(1.0659418) q[1];
sx q[1];
rz(-1.7168047) q[1];
sx q[1];
rz(2.8071383) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4518648) q[0];
sx q[0];
rz(-0.34669995) q[0];
sx q[0];
rz(2.6348216) q[0];
rz(0.43844079) q[2];
sx q[2];
rz(-1.4435569) q[2];
sx q[2];
rz(-1.3194989) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.85035861) q[1];
sx q[1];
rz(-2.6667074) q[1];
sx q[1];
rz(2.8104981) q[1];
rz(-pi) q[2];
rz(-0.45463698) q[3];
sx q[3];
rz(-2.2918923) q[3];
sx q[3];
rz(-2.7434512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7740384) q[2];
sx q[2];
rz(-0.84711051) q[2];
sx q[2];
rz(-0.96662194) q[2];
rz(-1.4878081) q[3];
sx q[3];
rz(-0.32854587) q[3];
sx q[3];
rz(3.0706792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8129355) q[0];
sx q[0];
rz(-0.85449496) q[0];
sx q[0];
rz(0.27035126) q[0];
rz(-1.5125754) q[1];
sx q[1];
rz(-0.83643475) q[1];
sx q[1];
rz(-0.62634748) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0369025) q[0];
sx q[0];
rz(-2.2459163) q[0];
sx q[0];
rz(2.7319607) q[0];
rz(-pi) q[1];
rz(0.43105189) q[2];
sx q[2];
rz(-0.97857514) q[2];
sx q[2];
rz(2.2811802) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8423564) q[1];
sx q[1];
rz(-1.5299712) q[1];
sx q[1];
rz(0.090260669) q[1];
rz(1.8966497) q[3];
sx q[3];
rz(-1.9385466) q[3];
sx q[3];
rz(2.8830584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7865929) q[2];
sx q[2];
rz(-2.2248416) q[2];
sx q[2];
rz(-1.5724486) q[2];
rz(-0.7507503) q[3];
sx q[3];
rz(-2.7995977) q[3];
sx q[3];
rz(-0.87740889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56562051) q[0];
sx q[0];
rz(-1.8436057) q[0];
sx q[0];
rz(2.5634503) q[0];
rz(-1.3427973) q[1];
sx q[1];
rz(-0.31675757) q[1];
sx q[1];
rz(-2.996179) q[1];
rz(1.5707468) q[2];
sx q[2];
rz(-1.4484828) q[2];
sx q[2];
rz(3.0640329) q[2];
rz(-0.20082898) q[3];
sx q[3];
rz(-0.26293892) q[3];
sx q[3];
rz(0.50504167) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
