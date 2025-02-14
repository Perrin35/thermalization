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
rz(-1.5768749) q[0];
sx q[0];
rz(-1.8160507) q[0];
sx q[0];
rz(-2.5817885) q[0];
rz(1.3979823) q[1];
sx q[1];
rz(-1.8140732) q[1];
sx q[1];
rz(-1.7550069) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5965393) q[0];
sx q[0];
rz(-2.5395416) q[0];
sx q[0];
rz(2.028378) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29560372) q[2];
sx q[2];
rz(-1.4628551) q[2];
sx q[2];
rz(2.2578277) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4475354) q[1];
sx q[1];
rz(-2.9172998) q[1];
sx q[1];
rz(-0.38557012) q[1];
rz(-pi) q[2];
rz(1.54744) q[3];
sx q[3];
rz(-0.29072116) q[3];
sx q[3];
rz(-0.61222044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5379415) q[2];
sx q[2];
rz(-2.5141022) q[2];
sx q[2];
rz(0.44169912) q[2];
rz(2.3407827) q[3];
sx q[3];
rz(-1.8947314) q[3];
sx q[3];
rz(2.7504564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(1.9073198) q[0];
sx q[0];
rz(-1.7664302) q[0];
sx q[0];
rz(-0.32670879) q[0];
rz(1.3455343) q[1];
sx q[1];
rz(-1.6163454) q[1];
sx q[1];
rz(-0.84652841) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.617793) q[0];
sx q[0];
rz(-1.1868068) q[0];
sx q[0];
rz(0.11597534) q[0];
rz(-0.88794215) q[2];
sx q[2];
rz(-1.3519716) q[2];
sx q[2];
rz(1.8162948) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.1607779) q[1];
sx q[1];
rz(-2.6123568) q[1];
sx q[1];
rz(1.5917626) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3896265) q[3];
sx q[3];
rz(-1.3526254) q[3];
sx q[3];
rz(-0.057540992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.232406) q[2];
sx q[2];
rz(-2.0441983) q[2];
sx q[2];
rz(-2.6573913) q[2];
rz(-1.8353315) q[3];
sx q[3];
rz(-0.56070119) q[3];
sx q[3];
rz(1.1854712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
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
rz(0.66121286) q[0];
sx q[0];
rz(-2.4112207) q[0];
sx q[0];
rz(-0.54006201) q[0];
rz(1.9909319) q[1];
sx q[1];
rz(-1.2723424) q[1];
sx q[1];
rz(-2.5475492) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68349379) q[0];
sx q[0];
rz(-0.20266315) q[0];
sx q[0];
rz(-2.1432102) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8833022) q[2];
sx q[2];
rz(-2.2021026) q[2];
sx q[2];
rz(0.7429276) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3004527) q[1];
sx q[1];
rz(-0.81592919) q[1];
sx q[1];
rz(-0.12187477) q[1];
x q[2];
rz(-2.0527127) q[3];
sx q[3];
rz(-1.6301117) q[3];
sx q[3];
rz(2.7869567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8192886) q[2];
sx q[2];
rz(-1.015859) q[2];
sx q[2];
rz(-2.4508396) q[2];
rz(-0.47762075) q[3];
sx q[3];
rz(-2.728929) q[3];
sx q[3];
rz(2.7121108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2531042) q[0];
sx q[0];
rz(-0.34585837) q[0];
sx q[0];
rz(0.75230569) q[0];
rz(-0.90657702) q[1];
sx q[1];
rz(-0.38175672) q[1];
sx q[1];
rz(-2.9125772) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54628263) q[0];
sx q[0];
rz(-0.37308274) q[0];
sx q[0];
rz(2.4855916) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3439517) q[2];
sx q[2];
rz(-1.0431492) q[2];
sx q[2];
rz(-2.2259852) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.467207) q[1];
sx q[1];
rz(-1.3027281) q[1];
sx q[1];
rz(1.3246791) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.205577) q[3];
sx q[3];
rz(-2.1250279) q[3];
sx q[3];
rz(-0.54868719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0166246) q[2];
sx q[2];
rz(-0.11398537) q[2];
sx q[2];
rz(2.1086741) q[2];
rz(-1.0659263) q[3];
sx q[3];
rz(-1.4706127) q[3];
sx q[3];
rz(1.8377931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2579987) q[0];
sx q[0];
rz(-0.086253919) q[0];
sx q[0];
rz(-1.8121383) q[0];
rz(0.52779245) q[1];
sx q[1];
rz(-1.3161491) q[1];
sx q[1];
rz(0.39906183) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51080158) q[0];
sx q[0];
rz(-0.36040655) q[0];
sx q[0];
rz(1.3494169) q[0];
x q[1];
rz(1.1588604) q[2];
sx q[2];
rz(-2.2592681) q[2];
sx q[2];
rz(2.7044123) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.1414772) q[1];
sx q[1];
rz(-2.2908604) q[1];
sx q[1];
rz(-1.9898371) q[1];
rz(-2.9513017) q[3];
sx q[3];
rz(-0.88493333) q[3];
sx q[3];
rz(-3.0113285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1479147) q[2];
sx q[2];
rz(-2.3892011) q[2];
sx q[2];
rz(-2.1368775) q[2];
rz(0.29804722) q[3];
sx q[3];
rz(-1.3532676) q[3];
sx q[3];
rz(-1.0684439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7816958) q[0];
sx q[0];
rz(-1.5157461) q[0];
sx q[0];
rz(-1.1915278) q[0];
rz(-0.53939348) q[1];
sx q[1];
rz(-1.0788147) q[1];
sx q[1];
rz(-2.3753812) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2711002) q[0];
sx q[0];
rz(-1.9166166) q[0];
sx q[0];
rz(2.150751) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4922392) q[2];
sx q[2];
rz(-2.1198065) q[2];
sx q[2];
rz(-1.8446814) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.898452) q[1];
sx q[1];
rz(-0.86075538) q[1];
sx q[1];
rz(2.8398127) q[1];
rz(-pi) q[2];
rz(-1.9412463) q[3];
sx q[3];
rz(-1.4806357) q[3];
sx q[3];
rz(-2.3925635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3300276) q[2];
sx q[2];
rz(-1.4469688) q[2];
sx q[2];
rz(2.2583029) q[2];
rz(-0.53089321) q[3];
sx q[3];
rz(-2.5575432) q[3];
sx q[3];
rz(-0.77597165) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74060488) q[0];
sx q[0];
rz(-0.29380909) q[0];
sx q[0];
rz(-0.64939943) q[0];
rz(-2.7679288) q[1];
sx q[1];
rz(-0.83761907) q[1];
sx q[1];
rz(-1.9992453) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23655791) q[0];
sx q[0];
rz(-1.8982186) q[0];
sx q[0];
rz(2.9740646) q[0];
rz(-pi) q[1];
rz(1.5863922) q[2];
sx q[2];
rz(-1.7806539) q[2];
sx q[2];
rz(-0.02053989) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.879952) q[1];
sx q[1];
rz(-1.8790207) q[1];
sx q[1];
rz(0.76517268) q[1];
x q[2];
rz(2.029923) q[3];
sx q[3];
rz(-1.6111501) q[3];
sx q[3];
rz(2.1488275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1335527) q[2];
sx q[2];
rz(-2.4456761) q[2];
sx q[2];
rz(2.7827061) q[2];
rz(-2.533203) q[3];
sx q[3];
rz(-1.5865934) q[3];
sx q[3];
rz(2.2801094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95978874) q[0];
sx q[0];
rz(-0.53878468) q[0];
sx q[0];
rz(-2.1754919) q[0];
rz(0.44202647) q[1];
sx q[1];
rz(-1.8531046) q[1];
sx q[1];
rz(-2.7332222) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.934487) q[0];
sx q[0];
rz(-1.4382595) q[0];
sx q[0];
rz(-0.026959628) q[0];
rz(-pi) q[1];
rz(2.8595905) q[2];
sx q[2];
rz(-0.15371727) q[2];
sx q[2];
rz(1.3957746) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.0023394452) q[1];
sx q[1];
rz(-2.1531257) q[1];
sx q[1];
rz(-2.9703388) q[1];
rz(-1.3134052) q[3];
sx q[3];
rz(-0.85581778) q[3];
sx q[3];
rz(-1.4090713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8926706) q[2];
sx q[2];
rz(-0.20888027) q[2];
sx q[2];
rz(-1.199031) q[2];
rz(1.8438953) q[3];
sx q[3];
rz(-1.0782995) q[3];
sx q[3];
rz(0.001999438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2099828) q[0];
sx q[0];
rz(-1.4298121) q[0];
sx q[0];
rz(2.0344875) q[0];
rz(-0.18868407) q[1];
sx q[1];
rz(-2.76177) q[1];
sx q[1];
rz(1.8170549) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4562603) q[0];
sx q[0];
rz(-0.59964824) q[0];
sx q[0];
rz(2.620282) q[0];
rz(1.5397838) q[2];
sx q[2];
rz(-2.5442269) q[2];
sx q[2];
rz(-1.7776521) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3834988) q[1];
sx q[1];
rz(-1.0304839) q[1];
sx q[1];
rz(-0.0030777085) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.56745133) q[3];
sx q[3];
rz(-2.3559915) q[3];
sx q[3];
rz(0.42450617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.3482427) q[2];
sx q[2];
rz(-0.57955727) q[2];
sx q[2];
rz(-2.9938475) q[2];
rz(-0.82792264) q[3];
sx q[3];
rz(-1.7235651) q[3];
sx q[3];
rz(-2.2016321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4350236) q[0];
sx q[0];
rz(-0.35025418) q[0];
sx q[0];
rz(1.5237923) q[0];
rz(-1.891547) q[1];
sx q[1];
rz(-2.3077272) q[1];
sx q[1];
rz(-2.7203383) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73250153) q[0];
sx q[0];
rz(-1.2048462) q[0];
sx q[0];
rz(0.41682314) q[0];
x q[1];
rz(-1.658861) q[2];
sx q[2];
rz(-0.95092623) q[2];
sx q[2];
rz(-1.5181092) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0388185) q[1];
sx q[1];
rz(-1.5613268) q[1];
sx q[1];
rz(0.086834309) q[1];
rz(-pi) q[2];
x q[2];
rz(0.47151819) q[3];
sx q[3];
rz(-2.6425411) q[3];
sx q[3];
rz(-1.3671444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1214352) q[2];
sx q[2];
rz(-2.406106) q[2];
sx q[2];
rz(2.6911733) q[2];
rz(-1.2823229) q[3];
sx q[3];
rz(-2.4560792) q[3];
sx q[3];
rz(-2.569017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0192075) q[0];
sx q[0];
rz(-1.5337802) q[0];
sx q[0];
rz(1.5564729) q[0];
rz(2.4284594) q[1];
sx q[1];
rz(-1.5819555) q[1];
sx q[1];
rz(-3.1176007) q[1];
rz(0.85999388) q[2];
sx q[2];
rz(-2.4888629) q[2];
sx q[2];
rz(0.98301907) q[2];
rz(-1.4982669) q[3];
sx q[3];
rz(-1.650643) q[3];
sx q[3];
rz(2.1908303) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
