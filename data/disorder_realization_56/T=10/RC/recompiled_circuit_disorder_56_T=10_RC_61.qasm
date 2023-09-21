OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.38874415) q[0];
sx q[0];
rz(3.677877) q[0];
sx q[0];
rz(10.372547) q[0];
rz(-1.3287969) q[1];
sx q[1];
rz(4.4089945) q[1];
sx q[1];
rz(10.452527) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7209137) q[0];
sx q[0];
rz(-1.5936216) q[0];
sx q[0];
rz(-2.7657397) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0601294) q[2];
sx q[2];
rz(-1.7667734) q[2];
sx q[2];
rz(0.34055647) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.23119588) q[1];
sx q[1];
rz(-0.80363552) q[1];
sx q[1];
rz(-2.5956144) q[1];
rz(-pi) q[2];
rz(-1.1231698) q[3];
sx q[3];
rz(-0.99222224) q[3];
sx q[3];
rz(2.6920126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0119005) q[2];
sx q[2];
rz(-1.7069858) q[2];
sx q[2];
rz(-1.0502846) q[2];
rz(1.1132647) q[3];
sx q[3];
rz(-0.89171019) q[3];
sx q[3];
rz(-3.0734857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.072409078) q[0];
sx q[0];
rz(-1.2453112) q[0];
sx q[0];
rz(0.29775277) q[0];
rz(2.521926) q[1];
sx q[1];
rz(-1.0071808) q[1];
sx q[1];
rz(2.0334977) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28716921) q[0];
sx q[0];
rz(-0.65432917) q[0];
sx q[0];
rz(1.9186583) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1531395) q[2];
sx q[2];
rz(-2.4465912) q[2];
sx q[2];
rz(-0.37441355) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.83924343) q[1];
sx q[1];
rz(-1.1855372) q[1];
sx q[1];
rz(-1.9531996) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8403345) q[3];
sx q[3];
rz(-0.88629913) q[3];
sx q[3];
rz(1.0052296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6790598) q[2];
sx q[2];
rz(-2.1647537) q[2];
sx q[2];
rz(2.1662946) q[2];
rz(2.2235928) q[3];
sx q[3];
rz(-1.8564329) q[3];
sx q[3];
rz(0.29618922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.179203) q[0];
sx q[0];
rz(-0.85150349) q[0];
sx q[0];
rz(-0.54291022) q[0];
rz(-2.2593598) q[1];
sx q[1];
rz(-2.0062607) q[1];
sx q[1];
rz(2.1767445) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0438) q[0];
sx q[0];
rz(-1.291073) q[0];
sx q[0];
rz(-2.8066737) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0028814) q[2];
sx q[2];
rz(-2.841946) q[2];
sx q[2];
rz(-1.9195997) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.112903) q[1];
sx q[1];
rz(-1.6351055) q[1];
sx q[1];
rz(1.4494004) q[1];
rz(-pi) q[2];
rz(1.1448907) q[3];
sx q[3];
rz(-0.45255462) q[3];
sx q[3];
rz(1.9354613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0470011) q[2];
sx q[2];
rz(-2.5307405) q[2];
sx q[2];
rz(1.1331406) q[2];
rz(-0.23162332) q[3];
sx q[3];
rz(-1.8685721) q[3];
sx q[3];
rz(2.384322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27947458) q[0];
sx q[0];
rz(-0.010443895) q[0];
sx q[0];
rz(1.7650771) q[0];
rz(0.51849413) q[1];
sx q[1];
rz(-1.8771749) q[1];
sx q[1];
rz(0.24212295) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0644558) q[0];
sx q[0];
rz(-2.4479439) q[0];
sx q[0];
rz(-2.8855188) q[0];
rz(-pi) q[1];
rz(0.02853407) q[2];
sx q[2];
rz(-2.7118073) q[2];
sx q[2];
rz(-2.1536749) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.081269216) q[1];
sx q[1];
rz(-0.4821061) q[1];
sx q[1];
rz(0.37270765) q[1];
rz(-2.7235892) q[3];
sx q[3];
rz(-1.6998708) q[3];
sx q[3];
rz(-1.9263182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.1233998) q[2];
sx q[2];
rz(-0.96863666) q[2];
sx q[2];
rz(-3.0467395) q[2];
rz(1.3421966) q[3];
sx q[3];
rz(-1.7443402) q[3];
sx q[3];
rz(-2.7024787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.78385329) q[0];
sx q[0];
rz(-1.3107603) q[0];
sx q[0];
rz(-0.36079303) q[0];
rz(-1.3882673) q[1];
sx q[1];
rz(-1.3307945) q[1];
sx q[1];
rz(1.1345908) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.40588) q[0];
sx q[0];
rz(-2.1942733) q[0];
sx q[0];
rz(2.1924125) q[0];
rz(2.2573932) q[2];
sx q[2];
rz(-2.6066385) q[2];
sx q[2];
rz(1.7810437) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7736588) q[1];
sx q[1];
rz(-1.7264257) q[1];
sx q[1];
rz(-0.44981287) q[1];
x q[2];
rz(-1.4868823) q[3];
sx q[3];
rz(-1.1353555) q[3];
sx q[3];
rz(-0.30121379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1363042) q[2];
sx q[2];
rz(-1.7207928) q[2];
sx q[2];
rz(-0.57265442) q[2];
rz(-0.92875656) q[3];
sx q[3];
rz(-0.52162617) q[3];
sx q[3];
rz(1.9740392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8144433) q[0];
sx q[0];
rz(-1.0961908) q[0];
sx q[0];
rz(1.8922528) q[0];
rz(1.2231187) q[1];
sx q[1];
rz(-1.5253116) q[1];
sx q[1];
rz(1.1522326) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6845247) q[0];
sx q[0];
rz(-0.10604924) q[0];
sx q[0];
rz(-1.4507136) q[0];
rz(-pi) q[1];
rz(-1.1587028) q[2];
sx q[2];
rz(-2.0372143) q[2];
sx q[2];
rz(-2.7127624) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7631543) q[1];
sx q[1];
rz(-1.3707146) q[1];
sx q[1];
rz(-1.4640019) q[1];
rz(-pi) q[2];
rz(-2.4310962) q[3];
sx q[3];
rz(-0.81516719) q[3];
sx q[3];
rz(1.1991771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.46889177) q[2];
sx q[2];
rz(-1.9273309) q[2];
sx q[2];
rz(-2.1506298) q[2];
rz(0.64783603) q[3];
sx q[3];
rz(-0.9655374) q[3];
sx q[3];
rz(-2.794054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8836442) q[0];
sx q[0];
rz(-2.914496) q[0];
sx q[0];
rz(-0.062285475) q[0];
rz(-2.9557872) q[1];
sx q[1];
rz(-1.4567016) q[1];
sx q[1];
rz(-2.7468162) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5236854) q[0];
sx q[0];
rz(-2.6751408) q[0];
sx q[0];
rz(2.594069) q[0];
x q[1];
rz(2.4418418) q[2];
sx q[2];
rz(-0.47669461) q[2];
sx q[2];
rz(-1.9112019) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.863184) q[1];
sx q[1];
rz(-1.8179968) q[1];
sx q[1];
rz(-0.43988887) q[1];
rz(2.9881334) q[3];
sx q[3];
rz(-1.5342661) q[3];
sx q[3];
rz(0.032144459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4930967) q[2];
sx q[2];
rz(-0.33005565) q[2];
sx q[2];
rz(0.27302343) q[2];
rz(-1.8388883) q[3];
sx q[3];
rz(-1.8282993) q[3];
sx q[3];
rz(2.8295512) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39772314) q[0];
sx q[0];
rz(-2.0070772) q[0];
sx q[0];
rz(-1.8564818) q[0];
rz(-1.5015191) q[1];
sx q[1];
rz(-1.7506426) q[1];
sx q[1];
rz(1.3407019) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96269755) q[0];
sx q[0];
rz(-2.6216051) q[0];
sx q[0];
rz(2.2682701) q[0];
rz(-pi) q[1];
rz(-2.4457744) q[2];
sx q[2];
rz(-1.6098445) q[2];
sx q[2];
rz(0.85862904) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9400085) q[1];
sx q[1];
rz(-1.219698) q[1];
sx q[1];
rz(1.1854978) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.31524661) q[3];
sx q[3];
rz(-1.9936221) q[3];
sx q[3];
rz(1.8002312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1061219) q[2];
sx q[2];
rz(-2.3351228) q[2];
sx q[2];
rz(1.0236053) q[2];
rz(-2.9566531) q[3];
sx q[3];
rz(-2.7513294) q[3];
sx q[3];
rz(0.24188724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0446562) q[0];
sx q[0];
rz(-2.1450295) q[0];
sx q[0];
rz(1.6212844) q[0];
rz(2.8114491) q[1];
sx q[1];
rz(-1.2076999) q[1];
sx q[1];
rz(0.83713371) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5569816) q[0];
sx q[0];
rz(-2.4298482) q[0];
sx q[0];
rz(2.8296489) q[0];
rz(-pi) q[1];
rz(-0.43948549) q[2];
sx q[2];
rz(-0.62289933) q[2];
sx q[2];
rz(-2.5650131) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.21048927) q[1];
sx q[1];
rz(-1.6753907) q[1];
sx q[1];
rz(-1.7235669) q[1];
rz(-pi) q[2];
rz(-2.8505441) q[3];
sx q[3];
rz(-1.0927199) q[3];
sx q[3];
rz(-0.92791286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.60124406) q[2];
sx q[2];
rz(-2.0566172) q[2];
sx q[2];
rz(2.9619651) q[2];
rz(-0.99572292) q[3];
sx q[3];
rz(-1.2487753) q[3];
sx q[3];
rz(1.3109591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-2.3778465) q[0];
sx q[0];
rz(-2.7959931) q[0];
sx q[0];
rz(-2.0843704) q[0];
rz(0.10748848) q[1];
sx q[1];
rz(-1.8881533) q[1];
sx q[1];
rz(-0.9799788) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36469034) q[0];
sx q[0];
rz(-1.5199465) q[0];
sx q[0];
rz(1.6965673) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0699632) q[2];
sx q[2];
rz(-2.0344779) q[2];
sx q[2];
rz(-2.0641363) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.500463) q[1];
sx q[1];
rz(-0.49912057) q[1];
sx q[1];
rz(-2.0874546) q[1];
rz(-pi) q[2];
rz(0.1006871) q[3];
sx q[3];
rz(-1.8929314) q[3];
sx q[3];
rz(1.5299357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8367299) q[2];
sx q[2];
rz(-1.9367846) q[2];
sx q[2];
rz(2.1137962) q[2];
rz(-1.3868388) q[3];
sx q[3];
rz(-1.8324865) q[3];
sx q[3];
rz(-0.28361472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86823157) q[0];
sx q[0];
rz(-1.0366476) q[0];
sx q[0];
rz(-1.114053) q[0];
rz(-0.83203075) q[1];
sx q[1];
rz(-2.6770626) q[1];
sx q[1];
rz(-2.4774036) q[1];
rz(-1.2912512) q[2];
sx q[2];
rz(-0.57689473) q[2];
sx q[2];
rz(2.9070791) q[2];
rz(2.5084393) q[3];
sx q[3];
rz(-2.5452151) q[3];
sx q[3];
rz(2.7635318) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];