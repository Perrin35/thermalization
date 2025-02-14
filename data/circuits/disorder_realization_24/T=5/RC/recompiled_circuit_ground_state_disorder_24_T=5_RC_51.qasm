OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.4234023) q[0];
sx q[0];
rz(-0.0027522491) q[0];
sx q[0];
rz(-1.0116853) q[0];
rz(-2.6818795) q[1];
sx q[1];
rz(-1.3016394) q[1];
sx q[1];
rz(0.17170061) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7778439) q[0];
sx q[0];
rz(-0.96723377) q[0];
sx q[0];
rz(1.4715172) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9471875) q[2];
sx q[2];
rz(-1.6792337) q[2];
sx q[2];
rz(2.6284503) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.57331177) q[1];
sx q[1];
rz(-0.75158316) q[1];
sx q[1];
rz(0.64935301) q[1];
x q[2];
rz(1.0461476) q[3];
sx q[3];
rz(-3.0339315) q[3];
sx q[3];
rz(-0.3357418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1823938) q[2];
sx q[2];
rz(-0.48754075) q[2];
sx q[2];
rz(1.7215151) q[2];
rz(-1.1848263) q[3];
sx q[3];
rz(-1.5337557) q[3];
sx q[3];
rz(-1.0678631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(3.0576393) q[0];
sx q[0];
rz(-1.6211442) q[0];
sx q[0];
rz(-0.49348304) q[0];
rz(2.3656942) q[1];
sx q[1];
rz(-2.6319365) q[1];
sx q[1];
rz(0.50481558) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0579266) q[0];
sx q[0];
rz(-2.2424881) q[0];
sx q[0];
rz(-2.2855285) q[0];
x q[1];
rz(1.8680598) q[2];
sx q[2];
rz(-2.0784162) q[2];
sx q[2];
rz(1.0647237) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.82693726) q[1];
sx q[1];
rz(-0.75410226) q[1];
sx q[1];
rz(-0.60390632) q[1];
rz(-0.22498954) q[3];
sx q[3];
rz(-0.27249042) q[3];
sx q[3];
rz(2.837473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8043171) q[2];
sx q[2];
rz(-1.8307468) q[2];
sx q[2];
rz(2.6709225) q[2];
rz(2.5110631) q[3];
sx q[3];
rz(-2.1157406) q[3];
sx q[3];
rz(1.7414909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0640963) q[0];
sx q[0];
rz(-3.016576) q[0];
sx q[0];
rz(0.65573829) q[0];
rz(-2.8302622) q[1];
sx q[1];
rz(-2.2475524) q[1];
sx q[1];
rz(-0.071455926) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8690259) q[0];
sx q[0];
rz(-1.2203958) q[0];
sx q[0];
rz(3.0883771) q[0];
rz(1.9391483) q[2];
sx q[2];
rz(-1.8873653) q[2];
sx q[2];
rz(-2.9376415) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8904222) q[1];
sx q[1];
rz(-2.648472) q[1];
sx q[1];
rz(-1.7095196) q[1];
rz(-pi) q[2];
rz(0.53160588) q[3];
sx q[3];
rz(-1.7544477) q[3];
sx q[3];
rz(-1.4830228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2283356) q[2];
sx q[2];
rz(-1.5568638) q[2];
sx q[2];
rz(-0.060001686) q[2];
rz(1.9289198) q[3];
sx q[3];
rz(-0.69827497) q[3];
sx q[3];
rz(-2.3771299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54670984) q[0];
sx q[0];
rz(-1.3285652) q[0];
sx q[0];
rz(-0.57719624) q[0];
rz(-0.82950854) q[1];
sx q[1];
rz(-2.0894158) q[1];
sx q[1];
rz(-0.83782354) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4234377) q[0];
sx q[0];
rz(-2.1912247) q[0];
sx q[0];
rz(-1.4761094) q[0];
x q[1];
rz(0.98060645) q[2];
sx q[2];
rz(-1.6211121) q[2];
sx q[2];
rz(2.3286164) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7260704) q[1];
sx q[1];
rz(-0.94968098) q[1];
sx q[1];
rz(0.46393053) q[1];
rz(-0.44562037) q[3];
sx q[3];
rz(-2.1975448) q[3];
sx q[3];
rz(-2.8901951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1159346) q[2];
sx q[2];
rz(-1.5237153) q[2];
sx q[2];
rz(1.9177829) q[2];
rz(-0.4661679) q[3];
sx q[3];
rz(-2.7397082) q[3];
sx q[3];
rz(0.60552067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8849628) q[0];
sx q[0];
rz(-1.4361359) q[0];
sx q[0];
rz(3.0404941) q[0];
rz(0.413232) q[1];
sx q[1];
rz(-2.7006472) q[1];
sx q[1];
rz(-2.0535927) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4306385) q[0];
sx q[0];
rz(-1.0885518) q[0];
sx q[0];
rz(-2.7706465) q[0];
x q[1];
rz(1.9740943) q[2];
sx q[2];
rz(-0.77770644) q[2];
sx q[2];
rz(2.4571927) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1697547) q[1];
sx q[1];
rz(-0.4167476) q[1];
sx q[1];
rz(-2.2107812) q[1];
x q[2];
rz(-2.4721778) q[3];
sx q[3];
rz(-1.6261887) q[3];
sx q[3];
rz(-0.15633571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2130412) q[2];
sx q[2];
rz(-2.2553359) q[2];
sx q[2];
rz(1.586033) q[2];
rz(-2.0689615) q[3];
sx q[3];
rz(-1.5242679) q[3];
sx q[3];
rz(3.0090289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17305408) q[0];
sx q[0];
rz(-2.2569077) q[0];
sx q[0];
rz(1.435085) q[0];
rz(-2.3834719) q[1];
sx q[1];
rz(-2.4393612) q[1];
sx q[1];
rz(-0.10163669) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2588239) q[0];
sx q[0];
rz(-0.92494915) q[0];
sx q[0];
rz(-0.96667883) q[0];
rz(-pi) q[1];
rz(-2.2287057) q[2];
sx q[2];
rz(-0.55333558) q[2];
sx q[2];
rz(2.0408415) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.96856252) q[1];
sx q[1];
rz(-2.8391339) q[1];
sx q[1];
rz(-0.36424251) q[1];
x q[2];
rz(-0.93433617) q[3];
sx q[3];
rz(-2.015967) q[3];
sx q[3];
rz(1.0220774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.3397843) q[2];
sx q[2];
rz(-1.1016176) q[2];
sx q[2];
rz(-1.3327117) q[2];
rz(2.0594275) q[3];
sx q[3];
rz(-0.75282955) q[3];
sx q[3];
rz(-1.4235206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69671714) q[0];
sx q[0];
rz(-1.2514021) q[0];
sx q[0];
rz(0.74309293) q[0];
rz(-0.7630868) q[1];
sx q[1];
rz(-0.63947314) q[1];
sx q[1];
rz(2.2230164) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0325286) q[0];
sx q[0];
rz(-1.5957513) q[0];
sx q[0];
rz(3.1335013) q[0];
rz(-pi) q[1];
rz(-0.95942504) q[2];
sx q[2];
rz(-2.583873) q[2];
sx q[2];
rz(0.096867933) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9715226) q[1];
sx q[1];
rz(-2.5405209) q[1];
sx q[1];
rz(-0.55723377) q[1];
rz(-pi) q[2];
x q[2];
rz(0.13217632) q[3];
sx q[3];
rz(-0.72849792) q[3];
sx q[3];
rz(0.26154172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.7243728) q[2];
sx q[2];
rz(-1.959274) q[2];
sx q[2];
rz(-2.7057538) q[2];
rz(-0.11442746) q[3];
sx q[3];
rz(-0.2468214) q[3];
sx q[3];
rz(2.54336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8082387) q[0];
sx q[0];
rz(-0.55753189) q[0];
sx q[0];
rz(-0.58407855) q[0];
rz(0.43560478) q[1];
sx q[1];
rz(-1.4608811) q[1];
sx q[1];
rz(1.6216507) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0348957) q[0];
sx q[0];
rz(-0.55904065) q[0];
sx q[0];
rz(-0.99026545) q[0];
rz(-1.1618091) q[2];
sx q[2];
rz(-1.3700546) q[2];
sx q[2];
rz(0.44524064) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4055109) q[1];
sx q[1];
rz(-2.9353432) q[1];
sx q[1];
rz(-3.095605) q[1];
rz(2.5831843) q[3];
sx q[3];
rz(-0.30277751) q[3];
sx q[3];
rz(1.6807792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.074177563) q[2];
sx q[2];
rz(-1.7895074) q[2];
sx q[2];
rz(-2.0254859) q[2];
rz(-1.762278) q[3];
sx q[3];
rz(-2.0383056) q[3];
sx q[3];
rz(-3.0232159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8922358) q[0];
sx q[0];
rz(-3.0638969) q[0];
sx q[0];
rz(2.351601) q[0];
rz(0.063591592) q[1];
sx q[1];
rz(-2.5345232) q[1];
sx q[1];
rz(2.7008609) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4249464) q[0];
sx q[0];
rz(-1.8688723) q[0];
sx q[0];
rz(1.2407202) q[0];
rz(-pi) q[1];
x q[1];
rz(0.67040261) q[2];
sx q[2];
rz(-0.60302654) q[2];
sx q[2];
rz(-0.5926026) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1524618) q[1];
sx q[1];
rz(-0.50639443) q[1];
sx q[1];
rz(1.2269065) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1001631) q[3];
sx q[3];
rz(-1.3815615) q[3];
sx q[3];
rz(1.3845058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8792087) q[2];
sx q[2];
rz(-0.79664207) q[2];
sx q[2];
rz(-2.4764496) q[2];
rz(0.77196676) q[3];
sx q[3];
rz(-0.77163458) q[3];
sx q[3];
rz(1.5099248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36122286) q[0];
sx q[0];
rz(-0.64798111) q[0];
sx q[0];
rz(1.004647) q[0];
rz(1.7338344) q[1];
sx q[1];
rz(-14/(3*pi)) q[1];
sx q[1];
rz(-1.2900603) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72899517) q[0];
sx q[0];
rz(-0.83461232) q[0];
sx q[0];
rz(-2.8060629) q[0];
rz(2.1707613) q[2];
sx q[2];
rz(-1.4024251) q[2];
sx q[2];
rz(1.6554993) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5651144) q[1];
sx q[1];
rz(-1.8449515) q[1];
sx q[1];
rz(1.3512011) q[1];
rz(1.8220002) q[3];
sx q[3];
rz(-0.61579865) q[3];
sx q[3];
rz(1.766966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1349692) q[2];
sx q[2];
rz(-1.1639736) q[2];
sx q[2];
rz(2.8912365) q[2];
rz(0.97744673) q[3];
sx q[3];
rz(-1.9340632) q[3];
sx q[3];
rz(1.4704963) q[3];
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
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95494315) q[0];
sx q[0];
rz(-1.432812) q[0];
sx q[0];
rz(-1.1695255) q[0];
rz(2.3969338) q[1];
sx q[1];
rz(-2.3458993) q[1];
sx q[1];
rz(0.69534272) q[1];
rz(-2.9982243) q[2];
sx q[2];
rz(-1.558254) q[2];
sx q[2];
rz(-2.8009453) q[2];
rz(1.0741735) q[3];
sx q[3];
rz(-1.1840829) q[3];
sx q[3];
rz(0.89331762) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
