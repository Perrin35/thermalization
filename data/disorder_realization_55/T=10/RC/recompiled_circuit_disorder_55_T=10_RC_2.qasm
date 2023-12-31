OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5947333) q[0];
sx q[0];
rz(-1.5164627) q[0];
sx q[0];
rz(-2.8773142) q[0];
rz(2.1677986) q[1];
sx q[1];
rz(-1.9314613) q[1];
sx q[1];
rz(2.4063453) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7289294) q[0];
sx q[0];
rz(-0.67617765) q[0];
sx q[0];
rz(-0.23763188) q[0];
x q[1];
rz(-2.5392883) q[2];
sx q[2];
rz(-0.77568433) q[2];
sx q[2];
rz(-1.3210981) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.6890251) q[1];
sx q[1];
rz(-1.8568294) q[1];
sx q[1];
rz(2.9806115) q[1];
rz(1.5988889) q[3];
sx q[3];
rz(-2.1510604) q[3];
sx q[3];
rz(-0.42682808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4720817) q[2];
sx q[2];
rz(-1.8410204) q[2];
sx q[2];
rz(1.1038587) q[2];
rz(-1.8707229) q[3];
sx q[3];
rz(-1.2277675) q[3];
sx q[3];
rz(-0.27403533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1141777) q[0];
sx q[0];
rz(-0.52863055) q[0];
sx q[0];
rz(-0.43637481) q[0];
rz(-2.6787058) q[1];
sx q[1];
rz(-1.0375689) q[1];
sx q[1];
rz(2.8754821) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5824495) q[0];
sx q[0];
rz(-0.83280116) q[0];
sx q[0];
rz(2.0931431) q[0];
rz(-pi) q[1];
rz(-1.8196462) q[2];
sx q[2];
rz(-0.58983931) q[2];
sx q[2];
rz(-0.8283386) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.74782158) q[1];
sx q[1];
rz(-0.95612477) q[1];
sx q[1];
rz(-0.6786896) q[1];
rz(-pi) q[2];
rz(0.9640785) q[3];
sx q[3];
rz(-1.420141) q[3];
sx q[3];
rz(-1.7863303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6767072) q[2];
sx q[2];
rz(-1.2767982) q[2];
sx q[2];
rz(-0.51149386) q[2];
rz(0.809508) q[3];
sx q[3];
rz(-1.5313238) q[3];
sx q[3];
rz(2.8539343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70616102) q[0];
sx q[0];
rz(-1.5478739) q[0];
sx q[0];
rz(0.92873746) q[0];
rz(1.4061032) q[1];
sx q[1];
rz(-2.4415253) q[1];
sx q[1];
rz(1.3471289) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72235332) q[0];
sx q[0];
rz(-1.1142715) q[0];
sx q[0];
rz(-1.2786352) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8560266) q[2];
sx q[2];
rz(-1.7226379) q[2];
sx q[2];
rz(-2.1053932) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.320142) q[1];
sx q[1];
rz(-2.4648033) q[1];
sx q[1];
rz(2.0435145) q[1];
rz(2.2778355) q[3];
sx q[3];
rz(-1.45544) q[3];
sx q[3];
rz(2.399721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9946263) q[2];
sx q[2];
rz(-1.7549843) q[2];
sx q[2];
rz(1.6195126) q[2];
rz(-2.8772723) q[3];
sx q[3];
rz(-1.0364573) q[3];
sx q[3];
rz(-2.6456397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79214823) q[0];
sx q[0];
rz(-1.9932207) q[0];
sx q[0];
rz(-0.96570063) q[0];
rz(-0.72215885) q[1];
sx q[1];
rz(-1.637371) q[1];
sx q[1];
rz(2.5818363) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9352919) q[0];
sx q[0];
rz(-1.9959404) q[0];
sx q[0];
rz(-1.4817119) q[0];
rz(-pi) q[1];
rz(0.83300029) q[2];
sx q[2];
rz(-0.90655316) q[2];
sx q[2];
rz(-0.26495648) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.30448118) q[1];
sx q[1];
rz(-1.5251535) q[1];
sx q[1];
rz(-0.95506217) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9066554) q[3];
sx q[3];
rz(-0.80312356) q[3];
sx q[3];
rz(-1.1266176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.42797783) q[2];
sx q[2];
rz(-1.5087936) q[2];
sx q[2];
rz(1.7170061) q[2];
rz(-0.26040855) q[3];
sx q[3];
rz(-1.3924761) q[3];
sx q[3];
rz(2.7105455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.150862) q[0];
sx q[0];
rz(-1.4980415) q[0];
sx q[0];
rz(0.36636233) q[0];
rz(-1.5953966) q[1];
sx q[1];
rz(-0.55391824) q[1];
sx q[1];
rz(-0.34367925) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2491964) q[0];
sx q[0];
rz(-2.2049892) q[0];
sx q[0];
rz(-0.87974324) q[0];
rz(-0.86954388) q[2];
sx q[2];
rz(-2.0336667) q[2];
sx q[2];
rz(-0.16823828) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1829454) q[1];
sx q[1];
rz(-1.4970386) q[1];
sx q[1];
rz(1.0720836) q[1];
rz(-pi) q[2];
rz(-2.214659) q[3];
sx q[3];
rz(-1.5878131) q[3];
sx q[3];
rz(3.0683558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7498103) q[2];
sx q[2];
rz(-0.85931531) q[2];
sx q[2];
rz(-3.0878477) q[2];
rz(-1.404445) q[3];
sx q[3];
rz(-2.6440547) q[3];
sx q[3];
rz(0.29156175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6546201) q[0];
sx q[0];
rz(-2.4916861) q[0];
sx q[0];
rz(-2.3858331) q[0];
rz(0.02515633) q[1];
sx q[1];
rz(-0.92725602) q[1];
sx q[1];
rz(0.25973928) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.320967) q[0];
sx q[0];
rz(-0.48829406) q[0];
sx q[0];
rz(-2.2886306) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2345215) q[2];
sx q[2];
rz(-1.5521126) q[2];
sx q[2];
rz(-0.74355723) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0745084) q[1];
sx q[1];
rz(-1.2423406) q[1];
sx q[1];
rz(-2.902608) q[1];
rz(-pi) q[2];
rz(-3.1236157) q[3];
sx q[3];
rz(-2.020105) q[3];
sx q[3];
rz(0.22389212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.48866895) q[2];
sx q[2];
rz(-0.80604625) q[2];
sx q[2];
rz(3.0409813) q[2];
rz(-2.9595024) q[3];
sx q[3];
rz(-2.263335) q[3];
sx q[3];
rz(1.3476868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.1473734) q[0];
sx q[0];
rz(-1.7213151) q[0];
sx q[0];
rz(-2.8314262) q[0];
rz(-2.639333) q[1];
sx q[1];
rz(-2.6175833) q[1];
sx q[1];
rz(0.60595864) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3340671) q[0];
sx q[0];
rz(-1.345732) q[0];
sx q[0];
rz(-0.62691697) q[0];
rz(-2.546054) q[2];
sx q[2];
rz(-2.7936802) q[2];
sx q[2];
rz(-2.4728342) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.088674) q[1];
sx q[1];
rz(-2.5175736) q[1];
sx q[1];
rz(-0.96159972) q[1];
x q[2];
rz(-0.80375399) q[3];
sx q[3];
rz(-2.440212) q[3];
sx q[3];
rz(2.3826117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.24017748) q[2];
sx q[2];
rz(-1.9577953) q[2];
sx q[2];
rz(-2.288738) q[2];
rz(-1.7715706) q[3];
sx q[3];
rz(-1.6829237) q[3];
sx q[3];
rz(0.20496932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9119499) q[0];
sx q[0];
rz(-0.59597534) q[0];
sx q[0];
rz(-1.461331) q[0];
rz(1.4029067) q[1];
sx q[1];
rz(-0.97424126) q[1];
sx q[1];
rz(3.0775552) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1916699) q[0];
sx q[0];
rz(-0.64282286) q[0];
sx q[0];
rz(1.886801) q[0];
rz(-pi) q[1];
rz(-2.089558) q[2];
sx q[2];
rz(-2.5737737) q[2];
sx q[2];
rz(2.0702814) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1053972) q[1];
sx q[1];
rz(-0.35696778) q[1];
sx q[1];
rz(-3.0371975) q[1];
x q[2];
rz(-1.4488892) q[3];
sx q[3];
rz(-0.31414437) q[3];
sx q[3];
rz(-1.4172518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.091207592) q[2];
sx q[2];
rz(-0.63082266) q[2];
sx q[2];
rz(1.5861661) q[2];
rz(-0.88820109) q[3];
sx q[3];
rz(-1.9017838) q[3];
sx q[3];
rz(0.92938882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9973307) q[0];
sx q[0];
rz(-1.0634796) q[0];
sx q[0];
rz(-1.7096747) q[0];
rz(-0.56888467) q[1];
sx q[1];
rz(-2.6061997) q[1];
sx q[1];
rz(-1.127839) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9040684) q[0];
sx q[0];
rz(-1.5868109) q[0];
sx q[0];
rz(-2.9928656) q[0];
rz(-pi) q[1];
rz(0.53290368) q[2];
sx q[2];
rz(-1.0120631) q[2];
sx q[2];
rz(-1.2677873) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7334426) q[1];
sx q[1];
rz(-0.24194716) q[1];
sx q[1];
rz(-1.5223632) q[1];
rz(-pi) q[2];
x q[2];
rz(0.3663775) q[3];
sx q[3];
rz(-1.0421703) q[3];
sx q[3];
rz(2.5322994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.52788064) q[2];
sx q[2];
rz(-0.85313672) q[2];
sx q[2];
rz(2.9679427) q[2];
rz(0.33637834) q[3];
sx q[3];
rz(-1.222638) q[3];
sx q[3];
rz(2.7500847) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7062475) q[0];
sx q[0];
rz(-2.5543537) q[0];
sx q[0];
rz(1.6760814) q[0];
rz(0.82410518) q[1];
sx q[1];
rz(-1.5287639) q[1];
sx q[1];
rz(2.5691659) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7596282) q[0];
sx q[0];
rz(-2.1158764) q[0];
sx q[0];
rz(0.99960021) q[0];
rz(-0.27406759) q[2];
sx q[2];
rz(-2.0147418) q[2];
sx q[2];
rz(-2.519671) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6483557) q[1];
sx q[1];
rz(-1.8093458) q[1];
sx q[1];
rz(1.7931213) q[1];
x q[2];
rz(-1.8626067) q[3];
sx q[3];
rz(-2.2581873) q[3];
sx q[3];
rz(-2.92047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3616025) q[2];
sx q[2];
rz(-2.6670691) q[2];
sx q[2];
rz(-0.53722107) q[2];
rz(-2.0843263) q[3];
sx q[3];
rz(-0.89151645) q[3];
sx q[3];
rz(-0.69361544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.854241) q[0];
sx q[0];
rz(-1.1653405) q[0];
sx q[0];
rz(-1.5821138) q[0];
rz(-2.6782425) q[1];
sx q[1];
rz(-2.2644823) q[1];
sx q[1];
rz(1.5092441) q[1];
rz(-2.2255185) q[2];
sx q[2];
rz(-1.6486042) q[2];
sx q[2];
rz(-0.83124607) q[2];
rz(-2.0588856) q[3];
sx q[3];
rz(-1.6255717) q[3];
sx q[3];
rz(-2.0902904) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
