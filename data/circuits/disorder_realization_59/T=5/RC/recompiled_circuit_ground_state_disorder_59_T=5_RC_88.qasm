OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.6049603) q[0];
sx q[0];
rz(-2.9807615) q[0];
sx q[0];
rz(0.39128006) q[0];
rz(3.0985576) q[1];
sx q[1];
rz(-1.6208836) q[1];
sx q[1];
rz(0.33837858) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7926377) q[0];
sx q[0];
rz(-0.39039877) q[0];
sx q[0];
rz(-0.81116311) q[0];
x q[1];
rz(-0.050655575) q[2];
sx q[2];
rz(-1.7850375) q[2];
sx q[2];
rz(-0.07344499) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0300301) q[1];
sx q[1];
rz(-1.7610368) q[1];
sx q[1];
rz(-0.89221402) q[1];
rz(-pi) q[2];
rz(0.24177052) q[3];
sx q[3];
rz(-0.42754506) q[3];
sx q[3];
rz(-1.1113785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.36838621) q[2];
sx q[2];
rz(-2.1619004) q[2];
sx q[2];
rz(1.0580019) q[2];
rz(0.10231415) q[3];
sx q[3];
rz(-1.5240069) q[3];
sx q[3];
rz(1.7412831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31829396) q[0];
sx q[0];
rz(-0.75495356) q[0];
sx q[0];
rz(-0.21389432) q[0];
rz(1.8077883) q[1];
sx q[1];
rz(-2.6413481) q[1];
sx q[1];
rz(0.28876367) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5014415) q[0];
sx q[0];
rz(-2.1955865) q[0];
sx q[0];
rz(2.9286317) q[0];
rz(-2.3043121) q[2];
sx q[2];
rz(-2.0990666) q[2];
sx q[2];
rz(-2.3579896) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5369479) q[1];
sx q[1];
rz(-1.6004452) q[1];
sx q[1];
rz(1.349959) q[1];
x q[2];
rz(1.9085732) q[3];
sx q[3];
rz(-1.984341) q[3];
sx q[3];
rz(0.097584399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5220962) q[2];
sx q[2];
rz(-2.3213826) q[2];
sx q[2];
rz(1.8572469) q[2];
rz(-1.8340825) q[3];
sx q[3];
rz(-1.8239572) q[3];
sx q[3];
rz(0.6984843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73806015) q[0];
sx q[0];
rz(-2.0750676) q[0];
sx q[0];
rz(2.2962978) q[0];
rz(0.90557939) q[1];
sx q[1];
rz(-0.60494345) q[1];
sx q[1];
rz(1.9639429) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15525165) q[0];
sx q[0];
rz(-0.097106783) q[0];
sx q[0];
rz(0.20669051) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3272456) q[2];
sx q[2];
rz(-2.4422514) q[2];
sx q[2];
rz(1.4171079) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3556087) q[1];
sx q[1];
rz(-0.61265433) q[1];
sx q[1];
rz(0.36566465) q[1];
x q[2];
rz(1.7813453) q[3];
sx q[3];
rz(-2.4857134) q[3];
sx q[3];
rz(2.0149503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5251069) q[2];
sx q[2];
rz(-0.944204) q[2];
sx q[2];
rz(0.13060972) q[2];
rz(-1.3579926) q[3];
sx q[3];
rz(-2.1328752) q[3];
sx q[3];
rz(1.2880026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
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
rz(-2.4865049) q[0];
sx q[0];
rz(-0.26432744) q[0];
sx q[0];
rz(-0.41410145) q[0];
rz(-0.86530238) q[1];
sx q[1];
rz(-1.2369786) q[1];
sx q[1];
rz(-0.6032595) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8651486) q[0];
sx q[0];
rz(-1.8521327) q[0];
sx q[0];
rz(1.9758609) q[0];
x q[1];
rz(0.19152051) q[2];
sx q[2];
rz(-0.3580123) q[2];
sx q[2];
rz(0.70068089) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.016876) q[1];
sx q[1];
rz(-1.5278597) q[1];
sx q[1];
rz(0.34833585) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0217085) q[3];
sx q[3];
rz(-1.881424) q[3];
sx q[3];
rz(-1.436123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6174751) q[2];
sx q[2];
rz(-1.5323428) q[2];
sx q[2];
rz(-0.65940801) q[2];
rz(-0.13255969) q[3];
sx q[3];
rz(-2.5949251) q[3];
sx q[3];
rz(-0.70422712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34683126) q[0];
sx q[0];
rz(-2.1857388) q[0];
sx q[0];
rz(-2.8723248) q[0];
rz(1.6412093) q[1];
sx q[1];
rz(-1.8625448) q[1];
sx q[1];
rz(-1.0795116) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0844042) q[0];
sx q[0];
rz(-1.3916172) q[0];
sx q[0];
rz(-2.4898131) q[0];
x q[1];
rz(-2.4731878) q[2];
sx q[2];
rz(-0.89819095) q[2];
sx q[2];
rz(2.2318537) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1126571) q[1];
sx q[1];
rz(-1.9017692) q[1];
sx q[1];
rz(-1.0616567) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7825384) q[3];
sx q[3];
rz(-1.8820888) q[3];
sx q[3];
rz(2.9632729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.72011224) q[2];
sx q[2];
rz(-0.58028996) q[2];
sx q[2];
rz(-2.271358) q[2];
rz(-1.3358215) q[3];
sx q[3];
rz(-1.8060874) q[3];
sx q[3];
rz(0.67572063) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0077591) q[0];
sx q[0];
rz(-0.02820153) q[0];
sx q[0];
rz(-2.5303685) q[0];
rz(2.4080343) q[1];
sx q[1];
rz(-1.4708142) q[1];
sx q[1];
rz(-0.38331097) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4188273) q[0];
sx q[0];
rz(-1.8486029) q[0];
sx q[0];
rz(0.72577839) q[0];
rz(0.82742847) q[2];
sx q[2];
rz(-1.4216004) q[2];
sx q[2];
rz(-2.0713446) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.94180471) q[1];
sx q[1];
rz(-2.2082303) q[1];
sx q[1];
rz(-2.9822523) q[1];
x q[2];
rz(0.46054764) q[3];
sx q[3];
rz(-0.78118284) q[3];
sx q[3];
rz(2.6278842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7395301) q[2];
sx q[2];
rz(-2.5867088) q[2];
sx q[2];
rz(1.0129207) q[2];
rz(0.5976451) q[3];
sx q[3];
rz(-1.4089855) q[3];
sx q[3];
rz(0.92596084) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0585854) q[0];
sx q[0];
rz(-2.7124131) q[0];
sx q[0];
rz(3.10566) q[0];
rz(-1.2220471) q[1];
sx q[1];
rz(-2.4655894) q[1];
sx q[1];
rz(2.8935208) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1897514) q[0];
sx q[0];
rz(-1.3782256) q[0];
sx q[0];
rz(-1.309881) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4291968) q[2];
sx q[2];
rz(-1.0842825) q[2];
sx q[2];
rz(1.6213191) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9315459) q[1];
sx q[1];
rz(-1.7341053) q[1];
sx q[1];
rz(1.943066) q[1];
x q[2];
rz(0.10493829) q[3];
sx q[3];
rz(-1.5682568) q[3];
sx q[3];
rz(0.030935098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0566569) q[2];
sx q[2];
rz(-2.0818905) q[2];
sx q[2];
rz(2.6094931) q[2];
rz(0.45446864) q[3];
sx q[3];
rz(-1.5458509) q[3];
sx q[3];
rz(-1.8475629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71659511) q[0];
sx q[0];
rz(-2.0918562) q[0];
sx q[0];
rz(-0.24018921) q[0];
rz(-2.5364618) q[1];
sx q[1];
rz(-1.7938219) q[1];
sx q[1];
rz(0.64632195) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28531269) q[0];
sx q[0];
rz(-2.3512771) q[0];
sx q[0];
rz(1.1423443) q[0];
rz(-pi) q[1];
rz(1.6679627) q[2];
sx q[2];
rz(-0.36280131) q[2];
sx q[2];
rz(-1.3902612) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.58144795) q[1];
sx q[1];
rz(-1.7755076) q[1];
sx q[1];
rz(2.678656) q[1];
x q[2];
rz(1.8490995) q[3];
sx q[3];
rz(-2.4670546) q[3];
sx q[3];
rz(1.9703678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1120844) q[2];
sx q[2];
rz(-1.6479475) q[2];
sx q[2];
rz(0.081681577) q[2];
rz(-0.84780848) q[3];
sx q[3];
rz(-0.40926465) q[3];
sx q[3];
rz(-2.9186287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8001556) q[0];
sx q[0];
rz(-0.84469047) q[0];
sx q[0];
rz(0.89168125) q[0];
rz(-1.2314697) q[1];
sx q[1];
rz(-2.6530118) q[1];
sx q[1];
rz(-0.60985342) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9187136) q[0];
sx q[0];
rz(-2.0628669) q[0];
sx q[0];
rz(-2.6124873) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.78471176) q[2];
sx q[2];
rz(-1.5605188) q[2];
sx q[2];
rz(2.1267872) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.6285204) q[1];
sx q[1];
rz(-1.2781464) q[1];
sx q[1];
rz(0.77192941) q[1];
x q[2];
rz(1.5242349) q[3];
sx q[3];
rz(-2.7961429) q[3];
sx q[3];
rz(0.54333401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8651198) q[2];
sx q[2];
rz(-2.1742994) q[2];
sx q[2];
rz(-0.27061948) q[2];
rz(-2.596415) q[3];
sx q[3];
rz(-1.8137648) q[3];
sx q[3];
rz(2.8934208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64428627) q[0];
sx q[0];
rz(-1.3077284) q[0];
sx q[0];
rz(-0.18648952) q[0];
rz(-1.2093774) q[1];
sx q[1];
rz(-1.3554363) q[1];
sx q[1];
rz(3.1219416) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6015897) q[0];
sx q[0];
rz(-0.9603017) q[0];
sx q[0];
rz(2.6363346) q[0];
x q[1];
rz(-2.8803359) q[2];
sx q[2];
rz(-1.831033) q[2];
sx q[2];
rz(-0.90554905) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3557089) q[1];
sx q[1];
rz(-1.2194835) q[1];
sx q[1];
rz(1.7339049) q[1];
rz(-0.65852209) q[3];
sx q[3];
rz(-0.88044518) q[3];
sx q[3];
rz(3.0999121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7095715) q[2];
sx q[2];
rz(-2.186128) q[2];
sx q[2];
rz(-1.5891937) q[2];
rz(3.0324557) q[3];
sx q[3];
rz(-1.2723294) q[3];
sx q[3];
rz(-0.2763589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5149665) q[0];
sx q[0];
rz(-1.4780541) q[0];
sx q[0];
rz(1.9274101) q[0];
rz(1.7264438) q[1];
sx q[1];
rz(-1.0217923) q[1];
sx q[1];
rz(1.6820977) q[1];
rz(1.8564275) q[2];
sx q[2];
rz(-1.141333) q[2];
sx q[2];
rz(-0.0091576413) q[2];
rz(-1.6676025) q[3];
sx q[3];
rz(-1.9514241) q[3];
sx q[3];
rz(-1.5497691) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
