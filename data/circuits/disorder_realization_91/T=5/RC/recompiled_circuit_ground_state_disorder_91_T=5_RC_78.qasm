OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.087091669) q[0];
sx q[0];
rz(-2.2890685) q[0];
sx q[0];
rz(-0.26242119) q[0];
rz(-0.30811319) q[1];
sx q[1];
rz(-0.83642712) q[1];
sx q[1];
rz(-3.1321373) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1149166) q[0];
sx q[0];
rz(-2.466188) q[0];
sx q[0];
rz(-2.3079655) q[0];
rz(-pi) q[1];
rz(0.043536206) q[2];
sx q[2];
rz(-1.6400717) q[2];
sx q[2];
rz(2.2486698) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.63543273) q[1];
sx q[1];
rz(-2.3316158) q[1];
sx q[1];
rz(2.2732938) q[1];
x q[2];
rz(-1.057714) q[3];
sx q[3];
rz(-2.3145967) q[3];
sx q[3];
rz(1.3199453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0414163) q[2];
sx q[2];
rz(-2.2569423) q[2];
sx q[2];
rz(2.6803988) q[2];
rz(-0.50208107) q[3];
sx q[3];
rz(-0.82379782) q[3];
sx q[3];
rz(-1.7320777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0111888) q[0];
sx q[0];
rz(-1.6809502) q[0];
sx q[0];
rz(-0.23656626) q[0];
rz(-0.81831167) q[1];
sx q[1];
rz(-1.2032443) q[1];
sx q[1];
rz(2.1990105) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0088975178) q[0];
sx q[0];
rz(-0.030936154) q[0];
sx q[0];
rz(-1.6583468) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0636395) q[2];
sx q[2];
rz(-0.13826577) q[2];
sx q[2];
rz(-2.0035715) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2494804) q[1];
sx q[1];
rz(-1.5693519) q[1];
sx q[1];
rz(2.4066928) q[1];
rz(-pi) q[2];
rz(2.0560741) q[3];
sx q[3];
rz(-1.266978) q[3];
sx q[3];
rz(1.7485545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4054823) q[2];
sx q[2];
rz(-0.63850275) q[2];
sx q[2];
rz(-2.45641) q[2];
rz(-2.3403366) q[3];
sx q[3];
rz(-1.5056491) q[3];
sx q[3];
rz(1.2509468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7291229) q[0];
sx q[0];
rz(-2.6710489) q[0];
sx q[0];
rz(2.1614918) q[0];
rz(0.12067548) q[1];
sx q[1];
rz(-1.1680892) q[1];
sx q[1];
rz(1.4792222) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6355094) q[0];
sx q[0];
rz(-1.732462) q[0];
sx q[0];
rz(-0.76753214) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.50522106) q[2];
sx q[2];
rz(-1.1596934) q[2];
sx q[2];
rz(-1.1875325) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.90872902) q[1];
sx q[1];
rz(-1.5759948) q[1];
sx q[1];
rz(0.0089110891) q[1];
rz(1.9068662) q[3];
sx q[3];
rz(-2.1343291) q[3];
sx q[3];
rz(2.270171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7145308) q[2];
sx q[2];
rz(-1.6468628) q[2];
sx q[2];
rz(2.96116) q[2];
rz(0.42029542) q[3];
sx q[3];
rz(-2.1642978) q[3];
sx q[3];
rz(0.54496566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9412398) q[0];
sx q[0];
rz(-1.350133) q[0];
sx q[0];
rz(3.0923162) q[0];
rz(2.409528) q[1];
sx q[1];
rz(-2.2923636) q[1];
sx q[1];
rz(1.7656309) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43983641) q[0];
sx q[0];
rz(-1.6876093) q[0];
sx q[0];
rz(2.9573836) q[0];
rz(-pi) q[1];
rz(2.92166) q[2];
sx q[2];
rz(-2.3637407) q[2];
sx q[2];
rz(-2.5663515) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.56616773) q[1];
sx q[1];
rz(-1.638104) q[1];
sx q[1];
rz(-2.7507902) q[1];
x q[2];
rz(2.2985372) q[3];
sx q[3];
rz(-0.84765654) q[3];
sx q[3];
rz(1.3031194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0338318) q[2];
sx q[2];
rz(-1.7511448) q[2];
sx q[2];
rz(-2.4793009) q[2];
rz(-1.3307339) q[3];
sx q[3];
rz(-2.3190506) q[3];
sx q[3];
rz(-1.8432157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0996284) q[0];
sx q[0];
rz(-0.30655107) q[0];
sx q[0];
rz(2.1339259) q[0];
rz(0.10534605) q[1];
sx q[1];
rz(-2.4844929) q[1];
sx q[1];
rz(-2.9109921) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14734257) q[0];
sx q[0];
rz(-1.2201637) q[0];
sx q[0];
rz(-2.7870534) q[0];
x q[1];
rz(-1.9318387) q[2];
sx q[2];
rz(-1.3461543) q[2];
sx q[2];
rz(1.063907) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.28921738) q[1];
sx q[1];
rz(-1.4061478) q[1];
sx q[1];
rz(-1.6573805) q[1];
x q[2];
rz(-2.9419961) q[3];
sx q[3];
rz(-0.65465876) q[3];
sx q[3];
rz(0.22206355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0577804) q[2];
sx q[2];
rz(-0.96692204) q[2];
sx q[2];
rz(0.18873611) q[2];
rz(-1.905929) q[3];
sx q[3];
rz(-0.50714791) q[3];
sx q[3];
rz(1.88131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1756303) q[0];
sx q[0];
rz(-1.2164793) q[0];
sx q[0];
rz(0.29773444) q[0];
rz(0.072602428) q[1];
sx q[1];
rz(-0.68521348) q[1];
sx q[1];
rz(-2.4593478) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0370343) q[0];
sx q[0];
rz(-0.68892043) q[0];
sx q[0];
rz(-1.0997186) q[0];
x q[1];
rz(-0.14839006) q[2];
sx q[2];
rz(-2.1990657) q[2];
sx q[2];
rz(-2.4895957) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.161769) q[1];
sx q[1];
rz(-1.7186972) q[1];
sx q[1];
rz(-1.921341) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0892049) q[3];
sx q[3];
rz(-2.4757407) q[3];
sx q[3];
rz(2.7107554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5630774) q[2];
sx q[2];
rz(-0.35334057) q[2];
sx q[2];
rz(0.087470857) q[2];
rz(0.89720094) q[3];
sx q[3];
rz(-1.8950491) q[3];
sx q[3];
rz(2.7520666) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.099365756) q[0];
sx q[0];
rz(-2.0392188) q[0];
sx q[0];
rz(1.9788096) q[0];
rz(-0.21576628) q[1];
sx q[1];
rz(-0.81753221) q[1];
sx q[1];
rz(-0.93528265) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5007846) q[0];
sx q[0];
rz(-1.3324059) q[0];
sx q[0];
rz(2.3604908) q[0];
rz(-pi) q[1];
rz(-1.8290524) q[2];
sx q[2];
rz(-1.7134981) q[2];
sx q[2];
rz(-1.1772798) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.079283) q[1];
sx q[1];
rz(-2.3543962) q[1];
sx q[1];
rz(1.6823549) q[1];
rz(-pi) q[2];
rz(-2.3386674) q[3];
sx q[3];
rz(-2.7028137) q[3];
sx q[3];
rz(0.9641274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0420578) q[2];
sx q[2];
rz(-2.4877986) q[2];
sx q[2];
rz(0.81134861) q[2];
rz(0.30284303) q[3];
sx q[3];
rz(-1.9767913) q[3];
sx q[3];
rz(2.5109049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52043668) q[0];
sx q[0];
rz(-2.8148837) q[0];
sx q[0];
rz(2.4503571) q[0];
rz(0.65903819) q[1];
sx q[1];
rz(-1.5182511) q[1];
sx q[1];
rz(2.5504327) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1189445) q[0];
sx q[0];
rz(-1.582358) q[0];
sx q[0];
rz(1.287879) q[0];
rz(-1.0724655) q[2];
sx q[2];
rz(-1.8231824) q[2];
sx q[2];
rz(-1.4344189) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2745251) q[1];
sx q[1];
rz(-2.3127416) q[1];
sx q[1];
rz(-1.1222003) q[1];
rz(-pi) q[2];
rz(3.0847315) q[3];
sx q[3];
rz(-1.2788075) q[3];
sx q[3];
rz(-0.1869456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9968694) q[2];
sx q[2];
rz(-0.99088061) q[2];
sx q[2];
rz(1.6807751) q[2];
rz(0.75374323) q[3];
sx q[3];
rz(-1.4494579) q[3];
sx q[3];
rz(2.6678705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6853471) q[0];
sx q[0];
rz(-1.101475) q[0];
sx q[0];
rz(-0.22612485) q[0];
rz(-1.6926758) q[1];
sx q[1];
rz(-0.96335226) q[1];
sx q[1];
rz(-1.0015782) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9739363) q[0];
sx q[0];
rz(-1.7154124) q[0];
sx q[0];
rz(0.58982106) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6999845) q[2];
sx q[2];
rz(-0.749365) q[2];
sx q[2];
rz(0.53560773) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.68943095) q[1];
sx q[1];
rz(-2.9254483) q[1];
sx q[1];
rz(-0.5139022) q[1];
rz(-pi) q[2];
x q[2];
rz(0.68437741) q[3];
sx q[3];
rz(-1.3357786) q[3];
sx q[3];
rz(-2.7647465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.61233026) q[2];
sx q[2];
rz(-2.3713106) q[2];
sx q[2];
rz(0.15753499) q[2];
rz(-2.9869288) q[3];
sx q[3];
rz(-1.1636584) q[3];
sx q[3];
rz(-0.4304339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4100819) q[0];
sx q[0];
rz(-1.9022576) q[0];
sx q[0];
rz(-0.70247689) q[0];
rz(2.6055873) q[1];
sx q[1];
rz(-0.82967007) q[1];
sx q[1];
rz(1.3943256) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.047612543) q[0];
sx q[0];
rz(-1.5450006) q[0];
sx q[0];
rz(-1.734478) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6800266) q[2];
sx q[2];
rz(-1.6642206) q[2];
sx q[2];
rz(-1.0051073) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6977384) q[1];
sx q[1];
rz(-1.8540253) q[1];
sx q[1];
rz(2.1069788) q[1];
x q[2];
rz(3.1411885) q[3];
sx q[3];
rz(-2.7260216) q[3];
sx q[3];
rz(-3.0311751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1348306) q[2];
sx q[2];
rz(-2.4595021) q[2];
sx q[2];
rz(1.6949867) q[2];
rz(0.26099482) q[3];
sx q[3];
rz(-1.010453) q[3];
sx q[3];
rz(1.235777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21988729) q[0];
sx q[0];
rz(-0.6548665) q[0];
sx q[0];
rz(-1.0153216) q[0];
rz(1.075853) q[1];
sx q[1];
rz(-1.8668108) q[1];
sx q[1];
rz(1.6783953) q[1];
rz(2.867472) q[2];
sx q[2];
rz(-2.7514396) q[2];
sx q[2];
rz(-0.20382602) q[2];
rz(0.45810926) q[3];
sx q[3];
rz(-0.81732133) q[3];
sx q[3];
rz(2.2504239) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
