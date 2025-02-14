OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5724343) q[0];
sx q[0];
rz(-2.1433266) q[0];
sx q[0];
rz(-0.38842595) q[0];
rz(-1.2216964) q[1];
sx q[1];
rz(-2.1887527) q[1];
sx q[1];
rz(0.0739007) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.088053) q[0];
sx q[0];
rz(-0.5660935) q[0];
sx q[0];
rz(-2.7527656) q[0];
rz(-pi) q[1];
rz(2.2227395) q[2];
sx q[2];
rz(-1.5519189) q[2];
sx q[2];
rz(0.23052653) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.072920567) q[1];
sx q[1];
rz(-2.34413) q[1];
sx q[1];
rz(1.0571521) q[1];
rz(-1.8201572) q[3];
sx q[3];
rz(-1.3624745) q[3];
sx q[3];
rz(1.1492532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.6834324) q[2];
sx q[2];
rz(-1.4417803) q[2];
sx q[2];
rz(-1.3275576) q[2];
rz(-0.96539998) q[3];
sx q[3];
rz(-0.5894956) q[3];
sx q[3];
rz(2.086967) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4210159) q[0];
sx q[0];
rz(-3.133931) q[0];
sx q[0];
rz(1.3910008) q[0];
rz(-1.1945456) q[1];
sx q[1];
rz(-0.60940131) q[1];
sx q[1];
rz(0.70558277) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5783516) q[0];
sx q[0];
rz(-0.8734428) q[0];
sx q[0];
rz(-2.2388377) q[0];
rz(2.0031571) q[2];
sx q[2];
rz(-0.92022824) q[2];
sx q[2];
rz(0.43076602) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3485364) q[1];
sx q[1];
rz(-2.4677241) q[1];
sx q[1];
rz(-2.1516031) q[1];
rz(-pi) q[2];
x q[2];
rz(0.67882244) q[3];
sx q[3];
rz(-0.72094432) q[3];
sx q[3];
rz(-1.3852505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0671063) q[2];
sx q[2];
rz(-1.6697845) q[2];
sx q[2];
rz(-0.86522317) q[2];
rz(-1.5557965) q[3];
sx q[3];
rz(-2.9977048) q[3];
sx q[3];
rz(-1.5417064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9928352) q[0];
sx q[0];
rz(-0.8388297) q[0];
sx q[0];
rz(-1.6194153) q[0];
rz(-1.7332227) q[1];
sx q[1];
rz(-2.759628) q[1];
sx q[1];
rz(2.9011889) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6344389) q[0];
sx q[0];
rz(-0.85570645) q[0];
sx q[0];
rz(-1.2583744) q[0];
rz(-2.9536134) q[2];
sx q[2];
rz(-1.1795292) q[2];
sx q[2];
rz(-0.94598929) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.63096365) q[1];
sx q[1];
rz(-1.6121702) q[1];
sx q[1];
rz(0.10527412) q[1];
rz(-pi) q[2];
rz(-0.48146507) q[3];
sx q[3];
rz(-1.6908619) q[3];
sx q[3];
rz(2.1019746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0092885) q[2];
sx q[2];
rz(-1.665364) q[2];
sx q[2];
rz(1.08584) q[2];
rz(0.050431937) q[3];
sx q[3];
rz(-1.4112873) q[3];
sx q[3];
rz(2.1850695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10629912) q[0];
sx q[0];
rz(-1.6989919) q[0];
sx q[0];
rz(-0.77044368) q[0];
rz(2.2162407) q[1];
sx q[1];
rz(-1.6567566) q[1];
sx q[1];
rz(0.83522183) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2352426) q[0];
sx q[0];
rz(-0.49717227) q[0];
sx q[0];
rz(-3.1072561) q[0];
rz(-pi) q[1];
rz(0.3496895) q[2];
sx q[2];
rz(-0.2564132) q[2];
sx q[2];
rz(1.0832322) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.28047565) q[1];
sx q[1];
rz(-2.0973699) q[1];
sx q[1];
rz(-0.014650281) q[1];
x q[2];
rz(-1.6862482) q[3];
sx q[3];
rz(-1.770062) q[3];
sx q[3];
rz(-1.5448567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.21021065) q[2];
sx q[2];
rz(-1.8385734) q[2];
sx q[2];
rz(-2.4118928) q[2];
rz(0.99916712) q[3];
sx q[3];
rz(-1.8347284) q[3];
sx q[3];
rz(2.6877747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5900742) q[0];
sx q[0];
rz(-2.9158264) q[0];
sx q[0];
rz(2.1424868) q[0];
rz(-2.0935811) q[1];
sx q[1];
rz(-2.4297355) q[1];
sx q[1];
rz(0.5447095) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7180965) q[0];
sx q[0];
rz(-2.831651) q[0];
sx q[0];
rz(0.62356068) q[0];
rz(0.72417132) q[2];
sx q[2];
rz(-1.0909547) q[2];
sx q[2];
rz(-1.9373231) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6259768) q[1];
sx q[1];
rz(-2.5624609) q[1];
sx q[1];
rz(-0.34767751) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7361511) q[3];
sx q[3];
rz(-1.5675987) q[3];
sx q[3];
rz(2.1759667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.58392945) q[2];
sx q[2];
rz(-1.1489392) q[2];
sx q[2];
rz(-1.0168797) q[2];
rz(2.478638) q[3];
sx q[3];
rz(-0.69279492) q[3];
sx q[3];
rz(1.3346765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3589631) q[0];
sx q[0];
rz(-2.6435659) q[0];
sx q[0];
rz(1.9899415) q[0];
rz(-0.89371124) q[1];
sx q[1];
rz(-1.7620554) q[1];
sx q[1];
rz(3.025257) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6457845) q[0];
sx q[0];
rz(-0.48775717) q[0];
sx q[0];
rz(1.8013823) q[0];
rz(-1.9151808) q[2];
sx q[2];
rz(-1.55416) q[2];
sx q[2];
rz(0.045147506) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.92116881) q[1];
sx q[1];
rz(-1.0225778) q[1];
sx q[1];
rz(2.9698644) q[1];
rz(-pi) q[2];
rz(0.086236091) q[3];
sx q[3];
rz(-2.3914118) q[3];
sx q[3];
rz(-2.3762109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5748888) q[2];
sx q[2];
rz(-2.7723007) q[2];
sx q[2];
rz(1.0075547) q[2];
rz(1.9762074) q[3];
sx q[3];
rz(-2.8663965) q[3];
sx q[3];
rz(2.5113441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4646869) q[0];
sx q[0];
rz(-2.6755264) q[0];
sx q[0];
rz(0.16978547) q[0];
rz(-0.67974293) q[1];
sx q[1];
rz(-0.83234537) q[1];
sx q[1];
rz(0.92175093) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0089192275) q[0];
sx q[0];
rz(-2.8555495) q[0];
sx q[0];
rz(2.689365) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5160732) q[2];
sx q[2];
rz(-1.8773019) q[2];
sx q[2];
rz(-0.79105575) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.55619986) q[1];
sx q[1];
rz(-2.0631644) q[1];
sx q[1];
rz(0.3979759) q[1];
rz(-pi) q[2];
rz(0.27382752) q[3];
sx q[3];
rz(-1.3845155) q[3];
sx q[3];
rz(-2.8485988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.13806954) q[2];
sx q[2];
rz(-1.4813083) q[2];
sx q[2];
rz(1.8180234) q[2];
rz(2.0924163) q[3];
sx q[3];
rz(-1.1314355) q[3];
sx q[3];
rz(1.3808892) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84999371) q[0];
sx q[0];
rz(-1.1350564) q[0];
sx q[0];
rz(2.3681613) q[0];
rz(-0.4666346) q[1];
sx q[1];
rz(-2.0687053) q[1];
sx q[1];
rz(3.1007865) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3916676) q[0];
sx q[0];
rz(-3.0081869) q[0];
sx q[0];
rz(1.6682503) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.62760533) q[2];
sx q[2];
rz(-0.96663953) q[2];
sx q[2];
rz(1.0221572) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1010712) q[1];
sx q[1];
rz(-0.98993976) q[1];
sx q[1];
rz(3.0809666) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0061408) q[3];
sx q[3];
rz(-1.7164383) q[3];
sx q[3];
rz(0.70051769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9890954) q[2];
sx q[2];
rz(-2.5941807) q[2];
sx q[2];
rz(-0.034817783) q[2];
rz(3.1133437) q[3];
sx q[3];
rz(-2.0939128) q[3];
sx q[3];
rz(-0.69407535) q[3];
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
rz(1.6377653) q[0];
sx q[0];
rz(-2.4389508) q[0];
sx q[0];
rz(1.0850061) q[0];
rz(-1.7833692) q[1];
sx q[1];
rz(-0.63301507) q[1];
sx q[1];
rz(-2.0620652) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87303783) q[0];
sx q[0];
rz(-1.7760906) q[0];
sx q[0];
rz(2.1389066) q[0];
x q[1];
rz(0.28330477) q[2];
sx q[2];
rz(-2.676801) q[2];
sx q[2];
rz(1.2086442) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4145045) q[1];
sx q[1];
rz(-1.8782756) q[1];
sx q[1];
rz(-0.25924637) q[1];
rz(-2.1302159) q[3];
sx q[3];
rz(-2.3314948) q[3];
sx q[3];
rz(1.427785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2435771) q[2];
sx q[2];
rz(-1.5960627) q[2];
sx q[2];
rz(0.071694516) q[2];
rz(0.60595766) q[3];
sx q[3];
rz(-0.76064435) q[3];
sx q[3];
rz(3.0683556) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83409413) q[0];
sx q[0];
rz(-1.5102757) q[0];
sx q[0];
rz(2.639005) q[0];
rz(1.7723627) q[1];
sx q[1];
rz(-1.2255729) q[1];
sx q[1];
rz(0.1756846) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0775107) q[0];
sx q[0];
rz(-0.33987576) q[0];
sx q[0];
rz(1.8301635) q[0];
rz(-pi) q[1];
rz(2.2403342) q[2];
sx q[2];
rz(-2.2961535) q[2];
sx q[2];
rz(2.8536316) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.60424544) q[1];
sx q[1];
rz(-2.2425695) q[1];
sx q[1];
rz(2.2293381) q[1];
x q[2];
rz(-2.0786344) q[3];
sx q[3];
rz(-0.70379095) q[3];
sx q[3];
rz(-2.2624571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.48675576) q[2];
sx q[2];
rz(-1.2590057) q[2];
sx q[2];
rz(2.7969825) q[2];
rz(-1.9136072) q[3];
sx q[3];
rz(-2.4495008) q[3];
sx q[3];
rz(2.1740348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6268613) q[0];
sx q[0];
rz(-1.3118962) q[0];
sx q[0];
rz(0.94919039) q[0];
rz(1.3728036) q[1];
sx q[1];
rz(-2.1887442) q[1];
sx q[1];
rz(-1.8123117) q[1];
rz(1.8362888) q[2];
sx q[2];
rz(-2.2659779) q[2];
sx q[2];
rz(-0.1018079) q[2];
rz(2.9459841) q[3];
sx q[3];
rz(-0.8630639) q[3];
sx q[3];
rz(0.56474781) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
