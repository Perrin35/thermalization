OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.91575032) q[0];
sx q[0];
rz(3.1728035) q[0];
sx q[0];
rz(6.7682545) q[0];
rz(-2.354061) q[1];
sx q[1];
rz(-2.1252316) q[1];
sx q[1];
rz(-2.7273942) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.480455) q[0];
sx q[0];
rz(-1.2167131) q[0];
sx q[0];
rz(2.7829091) q[0];
rz(-pi) q[1];
rz(1.7180874) q[2];
sx q[2];
rz(-0.5746595) q[2];
sx q[2];
rz(-1.6342083) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.58798446) q[1];
sx q[1];
rz(-1.506664) q[1];
sx q[1];
rz(-1.9967805) q[1];
rz(-pi) q[2];
rz(-1.4685417) q[3];
sx q[3];
rz(-1.3073321) q[3];
sx q[3];
rz(-1.9139569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2177314) q[2];
sx q[2];
rz(-1.2552746) q[2];
sx q[2];
rz(3.1100173) q[2];
rz(1.2565553) q[3];
sx q[3];
rz(-2.6387408) q[3];
sx q[3];
rz(-2.6149635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8388222) q[0];
sx q[0];
rz(-1.6844203) q[0];
sx q[0];
rz(2.9717428) q[0];
rz(-2.4376712) q[1];
sx q[1];
rz(-1.0715276) q[1];
sx q[1];
rz(0.53952113) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8115494) q[0];
sx q[0];
rz(-2.0199611) q[0];
sx q[0];
rz(-0.051785843) q[0];
rz(-0.45692921) q[2];
sx q[2];
rz(-0.77605844) q[2];
sx q[2];
rz(1.4571783) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0454228) q[1];
sx q[1];
rz(-1.3816557) q[1];
sx q[1];
rz(-1.063709) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.55595056) q[3];
sx q[3];
rz(-0.9405989) q[3];
sx q[3];
rz(-0.48660183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4743621) q[2];
sx q[2];
rz(-1.2381866) q[2];
sx q[2];
rz(2.898522) q[2];
rz(-0.66611755) q[3];
sx q[3];
rz(-2.5770498) q[3];
sx q[3];
rz(1.8977785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77984017) q[0];
sx q[0];
rz(-3.0267974) q[0];
sx q[0];
rz(2.6932122) q[0];
rz(1.386863) q[1];
sx q[1];
rz(-1.9877537) q[1];
sx q[1];
rz(0.2562491) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.019779531) q[0];
sx q[0];
rz(-0.67987961) q[0];
sx q[0];
rz(-2.7133184) q[0];
rz(-0.82661144) q[2];
sx q[2];
rz(-0.6140784) q[2];
sx q[2];
rz(-2.5342864) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5445404) q[1];
sx q[1];
rz(-1.8982732) q[1];
sx q[1];
rz(2.271133) q[1];
x q[2];
rz(-2.933421) q[3];
sx q[3];
rz(-0.18067193) q[3];
sx q[3];
rz(0.65073035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3391352) q[2];
sx q[2];
rz(-0.70636237) q[2];
sx q[2];
rz(0.57470542) q[2];
rz(-1.7859219) q[3];
sx q[3];
rz(-1.971743) q[3];
sx q[3];
rz(1.908196) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.213585) q[0];
sx q[0];
rz(-1.7127697) q[0];
sx q[0];
rz(2.8821049) q[0];
rz(1.150594) q[1];
sx q[1];
rz(-1.78803) q[1];
sx q[1];
rz(-0.73192275) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61011945) q[0];
sx q[0];
rz(-2.0391132) q[0];
sx q[0];
rz(-2.5584695) q[0];
rz(-pi) q[1];
rz(3.0078997) q[2];
sx q[2];
rz(-0.72172726) q[2];
sx q[2];
rz(1.114811) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.54938984) q[1];
sx q[1];
rz(-1.376774) q[1];
sx q[1];
rz(2.8787896) q[1];
rz(-pi) q[2];
rz(-0.51283522) q[3];
sx q[3];
rz(-2.8898015) q[3];
sx q[3];
rz(0.73392111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4449473) q[2];
sx q[2];
rz(-1.8680633) q[2];
sx q[2];
rz(-2.990492) q[2];
rz(2.5949196) q[3];
sx q[3];
rz(-2.0918545) q[3];
sx q[3];
rz(-2.7643519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54995173) q[0];
sx q[0];
rz(-1.0629835) q[0];
sx q[0];
rz(1.2623825) q[0];
rz(-1.4683912) q[1];
sx q[1];
rz(-0.60931283) q[1];
sx q[1];
rz(-0.79777065) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7724458) q[0];
sx q[0];
rz(-1.4506842) q[0];
sx q[0];
rz(0.0025047501) q[0];
rz(-2.5227929) q[2];
sx q[2];
rz(-0.36703645) q[2];
sx q[2];
rz(0.4180846) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5835727) q[1];
sx q[1];
rz(-0.94545525) q[1];
sx q[1];
rz(-0.11548345) q[1];
rz(-pi) q[2];
rz(1.2410774) q[3];
sx q[3];
rz(-1.9922678) q[3];
sx q[3];
rz(-0.7082522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.34565869) q[2];
sx q[2];
rz(-0.63085932) q[2];
sx q[2];
rz(-0.30203715) q[2];
rz(-1.1473514) q[3];
sx q[3];
rz(-1.6796422) q[3];
sx q[3];
rz(0.54774493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7323332) q[0];
sx q[0];
rz(-3.0497666) q[0];
sx q[0];
rz(-1.9858032) q[0];
rz(2.0571158) q[1];
sx q[1];
rz(-0.98025727) q[1];
sx q[1];
rz(0.070080431) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1228186) q[0];
sx q[0];
rz(-0.2304603) q[0];
sx q[0];
rz(-2.3019058) q[0];
x q[1];
rz(-2.5885133) q[2];
sx q[2];
rz(-2.2773909) q[2];
sx q[2];
rz(-1.4097708) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.97092123) q[1];
sx q[1];
rz(-1.0075924) q[1];
sx q[1];
rz(-1.7621653) q[1];
x q[2];
rz(1.1280941) q[3];
sx q[3];
rz(-0.08249313) q[3];
sx q[3];
rz(-0.024162956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8391116) q[2];
sx q[2];
rz(-1.182686) q[2];
sx q[2];
rz(-0.62136674) q[2];
rz(1.7012043) q[3];
sx q[3];
rz(-2.6337603) q[3];
sx q[3];
rz(2.5682209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086534111) q[0];
sx q[0];
rz(-1.4165514) q[0];
sx q[0];
rz(-2.281718) q[0];
rz(1.2043918) q[1];
sx q[1];
rz(-0.87221968) q[1];
sx q[1];
rz(0.0079356114) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36695776) q[0];
sx q[0];
rz(-1.8263706) q[0];
sx q[0];
rz(-2.5711683) q[0];
rz(-pi) q[1];
rz(-0.43086149) q[2];
sx q[2];
rz(-0.92509809) q[2];
sx q[2];
rz(0.25260392) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4850033) q[1];
sx q[1];
rz(-0.90066972) q[1];
sx q[1];
rz(-2.2336002) q[1];
x q[2];
rz(-2.8273724) q[3];
sx q[3];
rz(-1.1952956) q[3];
sx q[3];
rz(2.7995031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.69592151) q[2];
sx q[2];
rz(-1.9182703) q[2];
sx q[2];
rz(2.725214) q[2];
rz(-1.773206) q[3];
sx q[3];
rz(-1.2973283) q[3];
sx q[3];
rz(2.2369475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1798582) q[0];
sx q[0];
rz(-2.8680153) q[0];
sx q[0];
rz(-2.7767048) q[0];
rz(-2.2015613) q[1];
sx q[1];
rz(-0.53661984) q[1];
sx q[1];
rz(-1.6392802) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2973328) q[0];
sx q[0];
rz(-1.5607921) q[0];
sx q[0];
rz(0.25015932) q[0];
rz(-pi) q[1];
rz(2.084923) q[2];
sx q[2];
rz(-2.3377315) q[2];
sx q[2];
rz(-0.67509292) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9598436) q[1];
sx q[1];
rz(-1.0105003) q[1];
sx q[1];
rz(-0.77397857) q[1];
x q[2];
rz(1.7765462) q[3];
sx q[3];
rz(-1.4806517) q[3];
sx q[3];
rz(-2.3004325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6442948) q[2];
sx q[2];
rz(-0.50540322) q[2];
sx q[2];
rz(-1.0650744) q[2];
rz(-2.8403357) q[3];
sx q[3];
rz(-1.732429) q[3];
sx q[3];
rz(1.5208972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.168468) q[0];
sx q[0];
rz(-0.081806101) q[0];
sx q[0];
rz(-2.705943) q[0];
rz(-1.3849974) q[1];
sx q[1];
rz(-0.47043097) q[1];
sx q[1];
rz(2.7246144) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1736261) q[0];
sx q[0];
rz(-0.91714232) q[0];
sx q[0];
rz(3.0941512) q[0];
rz(-pi) q[1];
rz(0.68848227) q[2];
sx q[2];
rz(-2.3496029) q[2];
sx q[2];
rz(0.85306963) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4638085) q[1];
sx q[1];
rz(-2.6208372) q[1];
sx q[1];
rz(-1.2893454) q[1];
rz(-pi) q[2];
rz(-0.74420332) q[3];
sx q[3];
rz(-2.9508698) q[3];
sx q[3];
rz(0.71803367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0372662) q[2];
sx q[2];
rz(-1.6486847) q[2];
sx q[2];
rz(2.9157675) q[2];
rz(-0.2078235) q[3];
sx q[3];
rz(-2.4184629) q[3];
sx q[3];
rz(-2.5549755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41480961) q[0];
sx q[0];
rz(-0.2668969) q[0];
sx q[0];
rz(1.6171932) q[0];
rz(0.95364755) q[1];
sx q[1];
rz(-1.2748268) q[1];
sx q[1];
rz(1.8189925) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37227092) q[0];
sx q[0];
rz(-1.1936545) q[0];
sx q[0];
rz(-0.13052127) q[0];
x q[1];
rz(-0.25090353) q[2];
sx q[2];
rz(-2.0076027) q[2];
sx q[2];
rz(1.3383588) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9433371) q[1];
sx q[1];
rz(-2.2554734) q[1];
sx q[1];
rz(-0.68590045) q[1];
rz(-pi) q[2];
rz(2.9880452) q[3];
sx q[3];
rz(-2.2205177) q[3];
sx q[3];
rz(-2.6447907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3502729) q[2];
sx q[2];
rz(-1.046317) q[2];
sx q[2];
rz(1.0160758) q[2];
rz(-1.2223876) q[3];
sx q[3];
rz(-2.9639769) q[3];
sx q[3];
rz(-0.55541742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14324698) q[0];
sx q[0];
rz(-0.9220985) q[0];
sx q[0];
rz(-2.0621598) q[0];
rz(1.7779508) q[1];
sx q[1];
rz(-1.2294055) q[1];
sx q[1];
rz(1.3235863) q[1];
rz(-1.580711) q[2];
sx q[2];
rz(-2.0799939) q[2];
sx q[2];
rz(-1.6891198) q[2];
rz(1.2033403) q[3];
sx q[3];
rz(-0.24693476) q[3];
sx q[3];
rz(0.91528391) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];