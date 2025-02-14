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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.190024) q[0];
sx q[0];
rz(-0.61066884) q[0];
sx q[0];
rz(0.14279731) q[0];
rz(-pi) q[1];
rz(-1.4602939) q[2];
sx q[2];
rz(-1.7640451) q[2];
sx q[2];
rz(1.0363486) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.2312647) q[1];
sx q[1];
rz(-0.99580169) q[1];
sx q[1];
rz(2.0851589) q[1];
rz(-pi) q[2];
rz(3.0875077) q[3];
sx q[3];
rz(-1.6639317) q[3];
sx q[3];
rz(-2.2786811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.95919886) q[2];
sx q[2];
rz(-0.48754075) q[2];
sx q[2];
rz(1.7215151) q[2];
rz(1.1848263) q[3];
sx q[3];
rz(-1.607837) q[3];
sx q[3];
rz(2.0737295) q[3];
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
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083953388) q[0];
sx q[0];
rz(-1.6211442) q[0];
sx q[0];
rz(0.49348304) q[0];
rz(2.3656942) q[1];
sx q[1];
rz(-0.50965613) q[1];
sx q[1];
rz(-0.50481558) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0579266) q[0];
sx q[0];
rz(-0.89910451) q[0];
sx q[0];
rz(-0.85606411) q[0];
x q[1];
rz(1.8680598) q[2];
sx q[2];
rz(-2.0784162) q[2];
sx q[2];
rz(1.0647237) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9318092) q[1];
sx q[1];
rz(-1.9701013) q[1];
sx q[1];
rz(2.4834391) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8756469) q[3];
sx q[3];
rz(-1.6308745) q[3];
sx q[3];
rz(1.4836277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.33727553) q[2];
sx q[2];
rz(-1.8307468) q[2];
sx q[2];
rz(-2.6709225) q[2];
rz(2.5110631) q[3];
sx q[3];
rz(-2.1157406) q[3];
sx q[3];
rz(-1.4001018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0640963) q[0];
sx q[0];
rz(-3.016576) q[0];
sx q[0];
rz(-0.65573829) q[0];
rz(-0.31133044) q[1];
sx q[1];
rz(-0.89404023) q[1];
sx q[1];
rz(3.0701367) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27256672) q[0];
sx q[0];
rz(-1.9211968) q[0];
sx q[0];
rz(-0.053215543) q[0];
rz(-1.9391483) q[2];
sx q[2];
rz(-1.2542274) q[2];
sx q[2];
rz(0.20395111) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.944338) q[1];
sx q[1];
rz(-1.6363012) q[1];
sx q[1];
rz(-2.0599026) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7903665) q[3];
sx q[3];
rz(-0.5595419) q[3];
sx q[3];
rz(-2.9282687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2283356) q[2];
sx q[2];
rz(-1.5847289) q[2];
sx q[2];
rz(0.060001686) q[2];
rz(-1.9289198) q[3];
sx q[3];
rz(-0.69827497) q[3];
sx q[3];
rz(-0.76446271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54670984) q[0];
sx q[0];
rz(-1.3285652) q[0];
sx q[0];
rz(2.5643964) q[0];
rz(2.3120841) q[1];
sx q[1];
rz(-2.0894158) q[1];
sx q[1];
rz(-0.83782354) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71815495) q[0];
sx q[0];
rz(-2.1912247) q[0];
sx q[0];
rz(-1.4761094) q[0];
rz(2.1609862) q[2];
sx q[2];
rz(-1.5204805) q[2];
sx q[2];
rz(2.3286164) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7260704) q[1];
sx q[1];
rz(-0.94968098) q[1];
sx q[1];
rz(2.6776621) q[1];
x q[2];
rz(-2.2470705) q[3];
sx q[3];
rz(-1.2141879) q[3];
sx q[3];
rz(-1.5925533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.025658) q[2];
sx q[2];
rz(-1.5237153) q[2];
sx q[2];
rz(-1.9177829) q[2];
rz(-0.4661679) q[3];
sx q[3];
rz(-2.7397082) q[3];
sx q[3];
rz(0.60552067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25662988) q[0];
sx q[0];
rz(-1.4361359) q[0];
sx q[0];
rz(-3.0404941) q[0];
rz(2.7283607) q[1];
sx q[1];
rz(-0.44094545) q[1];
sx q[1];
rz(1.0879999) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0383069) q[0];
sx q[0];
rz(-1.8977471) q[0];
sx q[0];
rz(-2.0825544) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.83482439) q[2];
sx q[2];
rz(-1.8497648) q[2];
sx q[2];
rz(-0.5912515) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.99876304) q[1];
sx q[1];
rz(-1.8149477) q[1];
sx q[1];
rz(1.2296089) q[1];
rz(-1.641388) q[3];
sx q[3];
rz(-0.90259457) q[3];
sx q[3];
rz(-1.4582423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2130412) q[2];
sx q[2];
rz(-0.88625675) q[2];
sx q[2];
rz(-1.5555596) q[2];
rz(-1.0726311) q[3];
sx q[3];
rz(-1.5242679) q[3];
sx q[3];
rz(-3.0090289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(0.17305408) q[0];
sx q[0];
rz(-0.88468495) q[0];
sx q[0];
rz(1.7065077) q[0];
rz(0.75812078) q[1];
sx q[1];
rz(-0.70223141) q[1];
sx q[1];
rz(-3.039956) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0291189) q[0];
sx q[0];
rz(-2.2879507) q[0];
sx q[0];
rz(0.6458592) q[0];
x q[1];
rz(-1.1161713) q[2];
sx q[2];
rz(-1.897942) q[2];
sx q[2];
rz(1.0516372) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5530049) q[1];
sx q[1];
rz(-1.2887452) q[1];
sx q[1];
rz(-1.6815008) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2072565) q[3];
sx q[3];
rz(-2.015967) q[3];
sx q[3];
rz(-2.1195153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8018084) q[2];
sx q[2];
rz(-2.0399751) q[2];
sx q[2];
rz(-1.8088809) q[2];
rz(-1.0821651) q[3];
sx q[3];
rz(-0.75282955) q[3];
sx q[3];
rz(1.7180721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4448755) q[0];
sx q[0];
rz(-1.8901905) q[0];
sx q[0];
rz(-0.74309293) q[0];
rz(-2.3785059) q[1];
sx q[1];
rz(-2.5021195) q[1];
sx q[1];
rz(2.2230164) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4615304) q[0];
sx q[0];
rz(-1.5627075) q[0];
sx q[0];
rz(-1.5458406) q[0];
rz(-0.95942504) q[2];
sx q[2];
rz(-0.55771962) q[2];
sx q[2];
rz(-0.096867933) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2661675) q[1];
sx q[1];
rz(-1.2670749) q[1];
sx q[1];
rz(2.6145302) q[1];
rz(-pi) q[2];
rz(-3.0094163) q[3];
sx q[3];
rz(-0.72849792) q[3];
sx q[3];
rz(-2.8800509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.7243728) q[2];
sx q[2];
rz(-1.959274) q[2];
sx q[2];
rz(0.43583885) q[2];
rz(-0.11442746) q[3];
sx q[3];
rz(-2.8947713) q[3];
sx q[3];
rz(0.59823263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8082387) q[0];
sx q[0];
rz(-2.5840608) q[0];
sx q[0];
rz(-2.5575141) q[0];
rz(-0.43560478) q[1];
sx q[1];
rz(-1.4608811) q[1];
sx q[1];
rz(-1.6216507) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7652033) q[0];
sx q[0];
rz(-2.0302773) q[0];
sx q[0];
rz(-2.8110519) q[0];
x q[1];
rz(-1.1618091) q[2];
sx q[2];
rz(-1.3700546) q[2];
sx q[2];
rz(0.44524064) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.78306373) q[1];
sx q[1];
rz(-1.7768246) q[1];
sx q[1];
rz(-1.5804144) q[1];
rz(2.8826113) q[3];
sx q[3];
rz(-1.7294438) q[3];
sx q[3];
rz(2.4939031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0674151) q[2];
sx q[2];
rz(-1.7895074) q[2];
sx q[2];
rz(-1.1161067) q[2];
rz(-1.3793147) q[3];
sx q[3];
rz(-2.0383056) q[3];
sx q[3];
rz(-0.11837676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24935687) q[0];
sx q[0];
rz(-3.0638969) q[0];
sx q[0];
rz(-0.78999162) q[0];
rz(3.0780011) q[1];
sx q[1];
rz(-0.60706943) q[1];
sx q[1];
rz(2.7008609) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9544308) q[0];
sx q[0];
rz(-1.8858028) q[0];
sx q[0];
rz(0.31401547) q[0];
rz(-pi) q[1];
x q[1];
rz(2.47119) q[2];
sx q[2];
rz(-2.5385661) q[2];
sx q[2];
rz(2.5489901) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9891309) q[1];
sx q[1];
rz(-2.6351982) q[1];
sx q[1];
rz(1.9146862) q[1];
x q[2];
rz(-1.7601898) q[3];
sx q[3];
rz(-1.6114858) q[3];
sx q[3];
rz(0.19408801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8792087) q[2];
sx q[2];
rz(-0.79664207) q[2];
sx q[2];
rz(-0.66514307) q[2];
rz(-0.77196676) q[3];
sx q[3];
rz(-0.77163458) q[3];
sx q[3];
rz(1.6316679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36122286) q[0];
sx q[0];
rz(-0.64798111) q[0];
sx q[0];
rz(-1.004647) q[0];
rz(1.7338344) q[1];
sx q[1];
rz(-14/(3*pi)) q[1];
sx q[1];
rz(-1.2900603) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72899517) q[0];
sx q[0];
rz(-0.83461232) q[0];
sx q[0];
rz(-0.33552977) q[0];
rz(-0.97083135) q[2];
sx q[2];
rz(-1.7391676) q[2];
sx q[2];
rz(-1.6554993) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0755578) q[1];
sx q[1];
rz(-1.3595288) q[1];
sx q[1];
rz(0.28055625) q[1];
x q[2];
rz(-0.1741039) q[3];
sx q[3];
rz(-0.977036) q[3];
sx q[3];
rz(1.4623778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0066234) q[2];
sx q[2];
rz(-1.1639736) q[2];
sx q[2];
rz(-0.25035614) q[2];
rz(0.97744673) q[3];
sx q[3];
rz(-1.9340632) q[3];
sx q[3];
rz(1.4704963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95494315) q[0];
sx q[0];
rz(-1.7087806) q[0];
sx q[0];
rz(1.9720672) q[0];
rz(2.3969338) q[1];
sx q[1];
rz(-2.3458993) q[1];
sx q[1];
rz(0.69534272) q[1];
rz(3.0540285) q[2];
sx q[2];
rz(-2.9976805) q[2];
sx q[2];
rz(1.9981071) q[2];
rz(-0.43375265) q[3];
sx q[3];
rz(-2.0278145) q[3];
sx q[3];
rz(-0.47587004) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
