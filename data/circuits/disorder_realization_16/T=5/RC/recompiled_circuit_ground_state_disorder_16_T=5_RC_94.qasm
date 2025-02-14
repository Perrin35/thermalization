OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.0082173) q[0];
sx q[0];
rz(-1.2674588) q[0];
sx q[0];
rz(-0.01292364) q[0];
rz(0.68459964) q[1];
sx q[1];
rz(3.9404865) q[1];
sx q[1];
rz(10.482492) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81257129) q[0];
sx q[0];
rz(-0.91513915) q[0];
sx q[0];
rz(-1.283487) q[0];
rz(-1.0619668) q[2];
sx q[2];
rz(-0.86039174) q[2];
sx q[2];
rz(1.8343385) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.94016852) q[1];
sx q[1];
rz(-0.88646171) q[1];
sx q[1];
rz(-1.0971054) q[1];
rz(1.4679883) q[3];
sx q[3];
rz(-0.3936201) q[3];
sx q[3];
rz(3.032544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.58618033) q[2];
sx q[2];
rz(-1.5427898) q[2];
sx q[2];
rz(0.093322873) q[2];
rz(3.1210476) q[3];
sx q[3];
rz(-13/(16*pi)) q[3];
sx q[3];
rz(-1.3625905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(0.071863197) q[0];
sx q[0];
rz(-1.3988031) q[0];
sx q[0];
rz(-0.82759696) q[0];
rz(0.96356511) q[1];
sx q[1];
rz(-1.5255442) q[1];
sx q[1];
rz(0.73659426) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2824055) q[0];
sx q[0];
rz(-1.6466337) q[0];
sx q[0];
rz(-1.8453159) q[0];
rz(-pi) q[1];
x q[1];
rz(0.97073769) q[2];
sx q[2];
rz(-2.063437) q[2];
sx q[2];
rz(-3.0464095) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5399038) q[1];
sx q[1];
rz(-1.772953) q[1];
sx q[1];
rz(-1.2630995) q[1];
rz(-pi) q[2];
rz(-0.16421825) q[3];
sx q[3];
rz(-3.0053557) q[3];
sx q[3];
rz(3.1395903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5654512) q[2];
sx q[2];
rz(-2.5665923) q[2];
sx q[2];
rz(0.74419332) q[2];
rz(-0.60892504) q[3];
sx q[3];
rz(-0.78151339) q[3];
sx q[3];
rz(-1.6412546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7805444) q[0];
sx q[0];
rz(-0.86549509) q[0];
sx q[0];
rz(1.3737099) q[0];
rz(-2.3420077) q[1];
sx q[1];
rz(-2.1345963) q[1];
sx q[1];
rz(1.132157) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6705841) q[0];
sx q[0];
rz(-2.9078404) q[0];
sx q[0];
rz(-1.9712377) q[0];
rz(0.27081174) q[2];
sx q[2];
rz(-3.0584444) q[2];
sx q[2];
rz(-1.6390334) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.76946041) q[1];
sx q[1];
rz(-2.4677708) q[1];
sx q[1];
rz(-0.12427434) q[1];
x q[2];
rz(2.532652) q[3];
sx q[3];
rz(-1.534933) q[3];
sx q[3];
rz(-1.7882529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3867897) q[2];
sx q[2];
rz(-0.58480442) q[2];
sx q[2];
rz(-1.7542138) q[2];
rz(0.34902188) q[3];
sx q[3];
rz(-1.4533706) q[3];
sx q[3];
rz(2.6121228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3997407) q[0];
sx q[0];
rz(-0.96804237) q[0];
sx q[0];
rz(2.7834748) q[0];
rz(-2.070836) q[1];
sx q[1];
rz(-0.64859575) q[1];
sx q[1];
rz(-1.7020285) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5191906) q[0];
sx q[0];
rz(-0.82525142) q[0];
sx q[0];
rz(0.36143266) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.87936833) q[2];
sx q[2];
rz(-1.1739397) q[2];
sx q[2];
rz(0.59890998) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4434699) q[1];
sx q[1];
rz(-2.1204397) q[1];
sx q[1];
rz(2.2509031) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.374745) q[3];
sx q[3];
rz(-1.7550751) q[3];
sx q[3];
rz(-2.977918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4293356) q[2];
sx q[2];
rz(-1.2306932) q[2];
sx q[2];
rz(-0.62475359) q[2];
rz(1.5077) q[3];
sx q[3];
rz(-0.75993901) q[3];
sx q[3];
rz(-1.1225351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.7655012) q[0];
sx q[0];
rz(-0.88075817) q[0];
sx q[0];
rz(1.9951903) q[0];
rz(1.7064077) q[1];
sx q[1];
rz(-0.51992661) q[1];
sx q[1];
rz(0.82040876) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2489481) q[0];
sx q[0];
rz(-1.2995509) q[0];
sx q[0];
rz(0.091288996) q[0];
rz(-pi) q[1];
rz(-1.3455639) q[2];
sx q[2];
rz(-0.67906717) q[2];
sx q[2];
rz(-2.1175543) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8294404) q[1];
sx q[1];
rz(-1.7746801) q[1];
sx q[1];
rz(1.7704493) q[1];
rz(1.6040398) q[3];
sx q[3];
rz(-0.98131949) q[3];
sx q[3];
rz(-1.6575898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4121805) q[2];
sx q[2];
rz(-2.086144) q[2];
sx q[2];
rz(0.5385651) q[2];
rz(-1.4696848) q[3];
sx q[3];
rz(-1.4476176) q[3];
sx q[3];
rz(0.058549747) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3013714) q[0];
sx q[0];
rz(-3.0026307) q[0];
sx q[0];
rz(3.050991) q[0];
rz(2.7719356) q[1];
sx q[1];
rz(-1.548111) q[1];
sx q[1];
rz(2.7634117) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9210631) q[0];
sx q[0];
rz(-2.605633) q[0];
sx q[0];
rz(0.45951636) q[0];
rz(1.6126851) q[2];
sx q[2];
rz(-1.48078) q[2];
sx q[2];
rz(-0.57638327) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0269811) q[1];
sx q[1];
rz(-1.0798732) q[1];
sx q[1];
rz(-1.4308962) q[1];
x q[2];
rz(0.079433283) q[3];
sx q[3];
rz(-1.7609247) q[3];
sx q[3];
rz(3.0832689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.52046627) q[2];
sx q[2];
rz(-1.6660606) q[2];
sx q[2];
rz(3.0401518) q[2];
rz(-1.2927879) q[3];
sx q[3];
rz(-0.83270508) q[3];
sx q[3];
rz(1.7267936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3066278) q[0];
sx q[0];
rz(-1.2766301) q[0];
sx q[0];
rz(-0.90245885) q[0];
rz(2.9023671) q[1];
sx q[1];
rz(-1.7996412) q[1];
sx q[1];
rz(0.36144027) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28736605) q[0];
sx q[0];
rz(-2.0660095) q[0];
sx q[0];
rz(-0.12676858) q[0];
rz(-pi) q[1];
rz(-2.7290384) q[2];
sx q[2];
rz(-3.0027886) q[2];
sx q[2];
rz(-0.96497646) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.70009184) q[1];
sx q[1];
rz(-2.6569286) q[1];
sx q[1];
rz(-1.086471) q[1];
x q[2];
rz(0.42487259) q[3];
sx q[3];
rz(-0.36965224) q[3];
sx q[3];
rz(1.2583789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.07301894) q[2];
sx q[2];
rz(-1.8014427) q[2];
sx q[2];
rz(2.2646591) q[2];
rz(0.98313156) q[3];
sx q[3];
rz(-1.5638331) q[3];
sx q[3];
rz(-0.94299281) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32901949) q[0];
sx q[0];
rz(-2.946377) q[0];
sx q[0];
rz(-2.3028497) q[0];
rz(0.70873952) q[1];
sx q[1];
rz(-2.510431) q[1];
sx q[1];
rz(0.90604025) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86676407) q[0];
sx q[0];
rz(-1.5497009) q[0];
sx q[0];
rz(2.4167908) q[0];
rz(-pi) q[1];
rz(0.21561138) q[2];
sx q[2];
rz(-1.2279358) q[2];
sx q[2];
rz(-2.947399) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.38323572) q[1];
sx q[1];
rz(-0.89665186) q[1];
sx q[1];
rz(0.70887776) q[1];
x q[2];
rz(1.2165478) q[3];
sx q[3];
rz(-0.63874001) q[3];
sx q[3];
rz(0.50001345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.16443843) q[2];
sx q[2];
rz(-2.437037) q[2];
sx q[2];
rz(1.1158811) q[2];
rz(-2.5130533) q[3];
sx q[3];
rz(-1.081859) q[3];
sx q[3];
rz(2.5270497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3299265) q[0];
sx q[0];
rz(-0.8198494) q[0];
sx q[0];
rz(0.16127583) q[0];
rz(0.8423841) q[1];
sx q[1];
rz(-1.4295652) q[1];
sx q[1];
rz(-2.9715723) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5710167) q[0];
sx q[0];
rz(-1.5525155) q[0];
sx q[0];
rz(1.5626426) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4030946) q[2];
sx q[2];
rz(-1.1360886) q[2];
sx q[2];
rz(-2.0185883) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5304351) q[1];
sx q[1];
rz(-0.47947219) q[1];
sx q[1];
rz(-1.1392639) q[1];
rz(-1.5000072) q[3];
sx q[3];
rz(-0.38372358) q[3];
sx q[3];
rz(-3.0629223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1539479) q[2];
sx q[2];
rz(-1.5445856) q[2];
sx q[2];
rz(1.6551931) q[2];
rz(-3.0814643) q[3];
sx q[3];
rz(-2.3190053) q[3];
sx q[3];
rz(-1.7000343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5104367) q[0];
sx q[0];
rz(-0.031351723) q[0];
sx q[0];
rz(1.7106868) q[0];
rz(2.778964) q[1];
sx q[1];
rz(-1.6094094) q[1];
sx q[1];
rz(-1.3154202) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5595488) q[0];
sx q[0];
rz(-1.2799731) q[0];
sx q[0];
rz(-1.2469588) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9489297) q[2];
sx q[2];
rz(-1.6844498) q[2];
sx q[2];
rz(-3.0494351) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3676995) q[1];
sx q[1];
rz(-2.4871792) q[1];
sx q[1];
rz(-2.6952637) q[1];
rz(-pi) q[2];
x q[2];
rz(0.86851991) q[3];
sx q[3];
rz(-1.2099301) q[3];
sx q[3];
rz(-1.3165733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0843087) q[2];
sx q[2];
rz(-1.5949275) q[2];
sx q[2];
rz(-2.9810737) q[2];
rz(0.48216835) q[3];
sx q[3];
rz(-0.67626685) q[3];
sx q[3];
rz(-0.67960656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.01874825) q[0];
sx q[0];
rz(-1.3958805) q[0];
sx q[0];
rz(-0.91176283) q[0];
rz(1.979076) q[1];
sx q[1];
rz(-1.5717506) q[1];
sx q[1];
rz(2.7350978) q[1];
rz(0.087333655) q[2];
sx q[2];
rz(-1.6651911) q[2];
sx q[2];
rz(-1.390425) q[2];
rz(0.7708664) q[3];
sx q[3];
rz(-1.9996179) q[3];
sx q[3];
rz(-0.17948532) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
