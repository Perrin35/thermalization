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
rz(-2.8780706) q[0];
sx q[0];
rz(-1.4891013) q[0];
sx q[0];
rz(0.6058713) q[0];
x q[1];
rz(2.9471875) q[2];
sx q[2];
rz(-1.6792337) q[2];
sx q[2];
rz(-0.51314236) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.2312647) q[1];
sx q[1];
rz(-0.99580169) q[1];
sx q[1];
rz(1.0564338) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0954451) q[3];
sx q[3];
rz(-0.10766115) q[3];
sx q[3];
rz(-0.3357418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.95919886) q[2];
sx q[2];
rz(-0.48754075) q[2];
sx q[2];
rz(-1.4200776) q[2];
rz(-1.1848263) q[3];
sx q[3];
rz(-1.5337557) q[3];
sx q[3];
rz(2.0737295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0576393) q[0];
sx q[0];
rz(-1.6211442) q[0];
sx q[0];
rz(-2.6481096) q[0];
rz(-2.3656942) q[1];
sx q[1];
rz(-2.6319365) q[1];
sx q[1];
rz(2.6367771) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0579266) q[0];
sx q[0];
rz(-2.2424881) q[0];
sx q[0];
rz(0.85606411) q[0];
rz(-1.8680598) q[2];
sx q[2];
rz(-1.0631764) q[2];
sx q[2];
rz(1.0647237) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.82693726) q[1];
sx q[1];
rz(-0.75410226) q[1];
sx q[1];
rz(-0.60390632) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9166031) q[3];
sx q[3];
rz(-2.8691022) q[3];
sx q[3];
rz(2.837473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8043171) q[2];
sx q[2];
rz(-1.3108459) q[2];
sx q[2];
rz(2.6709225) q[2];
rz(2.5110631) q[3];
sx q[3];
rz(-1.0258521) q[3];
sx q[3];
rz(-1.7414909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.077496342) q[0];
sx q[0];
rz(-3.016576) q[0];
sx q[0];
rz(-0.65573829) q[0];
rz(2.8302622) q[1];
sx q[1];
rz(-2.2475524) q[1];
sx q[1];
rz(0.071455926) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0229683) q[0];
sx q[0];
rz(-0.35425348) q[0];
sx q[0];
rz(-1.4262761) q[0];
x q[1];
rz(2.3089789) q[2];
sx q[2];
rz(-2.660661) q[2];
sx q[2];
rz(-2.0455751) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1972547) q[1];
sx q[1];
rz(-1.6363012) q[1];
sx q[1];
rz(1.08169) q[1];
rz(1.7830332) q[3];
sx q[3];
rz(-2.0925412) q[3];
sx q[3];
rz(2.9468342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.91325703) q[2];
sx q[2];
rz(-1.5568638) q[2];
sx q[2];
rz(-3.081591) q[2];
rz(1.2126728) q[3];
sx q[3];
rz(-0.69827497) q[3];
sx q[3];
rz(-0.76446271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54670984) q[0];
sx q[0];
rz(-1.3285652) q[0];
sx q[0];
rz(-2.5643964) q[0];
rz(2.3120841) q[1];
sx q[1];
rz(-1.0521768) q[1];
sx q[1];
rz(-2.3037691) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5853607) q[0];
sx q[0];
rz(-0.62667003) q[0];
sx q[0];
rz(-3.0100432) q[0];
rz(-1.4805541) q[2];
sx q[2];
rz(-0.59207661) q[2];
sx q[2];
rz(0.83276487) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.41552222) q[1];
sx q[1];
rz(-0.94968098) q[1];
sx q[1];
rz(-2.6776621) q[1];
rz(-pi) q[2];
rz(-2.2470705) q[3];
sx q[3];
rz(-1.2141879) q[3];
sx q[3];
rz(1.5490393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.025658) q[2];
sx q[2];
rz(-1.6178774) q[2];
sx q[2];
rz(1.2238097) q[2];
rz(2.6754248) q[3];
sx q[3];
rz(-2.7397082) q[3];
sx q[3];
rz(-2.536072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25662988) q[0];
sx q[0];
rz(-1.7054568) q[0];
sx q[0];
rz(0.10109854) q[0];
rz(-2.7283607) q[1];
sx q[1];
rz(-2.7006472) q[1];
sx q[1];
rz(-2.0535927) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4306385) q[0];
sx q[0];
rz(-1.0885518) q[0];
sx q[0];
rz(-0.37094613) q[0];
rz(-pi) q[1];
rz(0.36878196) q[2];
sx q[2];
rz(-0.86930828) q[2];
sx q[2];
rz(-1.9176287) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6551825) q[1];
sx q[1];
rz(-1.2401199) q[1];
sx q[1];
rz(0.25844708) q[1];
rz(-pi) q[2];
rz(-2.4721778) q[3];
sx q[3];
rz(-1.515404) q[3];
sx q[3];
rz(0.15633571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.92855144) q[2];
sx q[2];
rz(-0.88625675) q[2];
sx q[2];
rz(1.586033) q[2];
rz(2.0689615) q[3];
sx q[3];
rz(-1.5242679) q[3];
sx q[3];
rz(0.13256375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17305408) q[0];
sx q[0];
rz(-2.2569077) q[0];
sx q[0];
rz(-1.435085) q[0];
rz(-2.3834719) q[1];
sx q[1];
rz(-2.4393612) q[1];
sx q[1];
rz(3.039956) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91824965) q[0];
sx q[0];
rz(-2.0416284) q[0];
sx q[0];
rz(2.4001394) q[0];
rz(2.2287057) q[2];
sx q[2];
rz(-0.55333558) q[2];
sx q[2];
rz(-2.0408415) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5530049) q[1];
sx q[1];
rz(-1.8528474) q[1];
sx q[1];
rz(1.4600919) q[1];
x q[2];
rz(-2.2472081) q[3];
sx q[3];
rz(-2.3830288) q[3];
sx q[3];
rz(1.0762843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.3397843) q[2];
sx q[2];
rz(-1.1016176) q[2];
sx q[2];
rz(-1.3327117) q[2];
rz(1.0821651) q[3];
sx q[3];
rz(-2.3887631) q[3];
sx q[3];
rz(-1.4235206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4448755) q[0];
sx q[0];
rz(-1.2514021) q[0];
sx q[0];
rz(-0.74309293) q[0];
rz(0.7630868) q[1];
sx q[1];
rz(-0.63947314) q[1];
sx q[1];
rz(0.9185763) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10906405) q[0];
sx q[0];
rz(-1.5957513) q[0];
sx q[0];
rz(-0.0080913261) q[0];
rz(-pi) q[1];
rz(-2.7977711) q[2];
sx q[2];
rz(-2.0190416) q[2];
sx q[2];
rz(-0.78735414) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2661675) q[1];
sx q[1];
rz(-1.2670749) q[1];
sx q[1];
rz(-2.6145302) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0094163) q[3];
sx q[3];
rz(-0.72849792) q[3];
sx q[3];
rz(0.26154172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4172198) q[2];
sx q[2];
rz(-1.1823187) q[2];
sx q[2];
rz(-2.7057538) q[2];
rz(3.0271652) q[3];
sx q[3];
rz(-0.2468214) q[3];
sx q[3];
rz(2.54336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8082387) q[0];
sx q[0];
rz(-2.5840608) q[0];
sx q[0];
rz(2.5575141) q[0];
rz(-0.43560478) q[1];
sx q[1];
rz(-1.4608811) q[1];
sx q[1];
rz(-1.6216507) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95660644) q[0];
sx q[0];
rz(-1.865956) q[0];
sx q[0];
rz(-1.0888238) q[0];
rz(-pi) q[1];
rz(2.043739) q[2];
sx q[2];
rz(-0.45308896) q[2];
sx q[2];
rz(2.447213) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3518924) q[1];
sx q[1];
rz(-1.580211) q[1];
sx q[1];
rz(-0.20603754) q[1];
rz(-pi) q[2];
rz(0.55840839) q[3];
sx q[3];
rz(-0.30277751) q[3];
sx q[3];
rz(-1.6807792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.074177563) q[2];
sx q[2];
rz(-1.7895074) q[2];
sx q[2];
rz(2.0254859) q[2];
rz(-1.3793147) q[3];
sx q[3];
rz(-2.0383056) q[3];
sx q[3];
rz(-0.11837676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8922358) q[0];
sx q[0];
rz(-0.07769575) q[0];
sx q[0];
rz(-0.78999162) q[0];
rz(3.0780011) q[1];
sx q[1];
rz(-0.60706943) q[1];
sx q[1];
rz(-0.44073179) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9544308) q[0];
sx q[0];
rz(-1.2557898) q[0];
sx q[0];
rz(-0.31401547) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1665383) q[2];
sx q[2];
rz(-2.031293) q[2];
sx q[2];
rz(-1.782589) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9891309) q[1];
sx q[1];
rz(-0.50639443) q[1];
sx q[1];
rz(1.9146862) q[1];
rz(-pi) q[2];
rz(1.7837672) q[3];
sx q[3];
rz(-0.19366385) q[3];
sx q[3];
rz(1.97399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.26238394) q[2];
sx q[2];
rz(-0.79664207) q[2];
sx q[2];
rz(2.4764496) q[2];
rz(-0.77196676) q[3];
sx q[3];
rz(-0.77163458) q[3];
sx q[3];
rz(1.6316679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
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
rz(-1.4077582) q[1];
sx q[1];
rz(-1.6561457) q[1];
sx q[1];
rz(1.2900603) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8915931) q[0];
sx q[0];
rz(-2.3457922) q[0];
sx q[0];
rz(1.222247) q[0];
x q[1];
rz(2.1707613) q[2];
sx q[2];
rz(-1.4024251) q[2];
sx q[2];
rz(1.6554993) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5764783) q[1];
sx q[1];
rz(-1.2966411) q[1];
sx q[1];
rz(1.7903916) q[1];
rz(-1.8220002) q[3];
sx q[3];
rz(-0.61579865) q[3];
sx q[3];
rz(1.3746266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1349692) q[2];
sx q[2];
rz(-1.9776191) q[2];
sx q[2];
rz(-0.25035614) q[2];
rz(-0.97744673) q[3];
sx q[3];
rz(-1.2075295) q[3];
sx q[3];
rz(1.4704963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95494315) q[0];
sx q[0];
rz(-1.7087806) q[0];
sx q[0];
rz(1.9720672) q[0];
rz(-0.74465887) q[1];
sx q[1];
rz(-2.3458993) q[1];
sx q[1];
rz(0.69534272) q[1];
rz(-2.9982243) q[2];
sx q[2];
rz(-1.558254) q[2];
sx q[2];
rz(-2.8009453) q[2];
rz(2.2779989) q[3];
sx q[3];
rz(-2.5222688) q[3];
sx q[3];
rz(1.8561192) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
