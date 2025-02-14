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
rz(0.45971316) q[1];
sx q[1];
rz(4.4432321) q[1];
sx q[1];
rz(9.2530773) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8780706) q[0];
sx q[0];
rz(-1.6524914) q[0];
sx q[0];
rz(0.6058713) q[0];
rz(-pi) q[1];
x q[1];
rz(0.19440513) q[2];
sx q[2];
rz(-1.6792337) q[2];
sx q[2];
rz(0.51314236) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.2312647) q[1];
sx q[1];
rz(-0.99580169) q[1];
sx q[1];
rz(-2.0851589) q[1];
rz(3.0875077) q[3];
sx q[3];
rz(-1.4776609) q[3];
sx q[3];
rz(2.2786811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.95919886) q[2];
sx q[2];
rz(-0.48754075) q[2];
sx q[2];
rz(-1.4200776) q[2];
rz(-1.1848263) q[3];
sx q[3];
rz(-1.5337557) q[3];
sx q[3];
rz(-1.0678631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0576393) q[0];
sx q[0];
rz(-1.5204484) q[0];
sx q[0];
rz(2.6481096) q[0];
rz(-0.77589846) q[1];
sx q[1];
rz(-0.50965613) q[1];
sx q[1];
rz(2.6367771) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.017708721) q[0];
sx q[0];
rz(-2.1095182) q[0];
sx q[0];
rz(-0.81102844) q[0];
x q[1];
rz(-0.4846862) q[2];
sx q[2];
rz(-0.58161608) q[2];
sx q[2];
rz(-1.6270552) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.82693726) q[1];
sx q[1];
rz(-2.3874904) q[1];
sx q[1];
rz(0.60390632) q[1];
x q[2];
rz(-0.26594578) q[3];
sx q[3];
rz(-1.6308745) q[3];
sx q[3];
rz(-1.4836277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.33727553) q[2];
sx q[2];
rz(-1.8307468) q[2];
sx q[2];
rz(0.47067019) q[2];
rz(-2.5110631) q[3];
sx q[3];
rz(-2.1157406) q[3];
sx q[3];
rz(-1.7414909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
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
rz(-3.0701367) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1186244) q[0];
sx q[0];
rz(-2.7873392) q[0];
sx q[0];
rz(-1.7153165) q[0];
rz(-2.803903) q[2];
sx q[2];
rz(-1.2215541) q[2];
sx q[2];
rz(1.8943292) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.944338) q[1];
sx q[1];
rz(-1.5052915) q[1];
sx q[1];
rz(1.08169) q[1];
x q[2];
rz(-1.7830332) q[3];
sx q[3];
rz(-2.0925412) q[3];
sx q[3];
rz(-2.9468342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2283356) q[2];
sx q[2];
rz(-1.5847289) q[2];
sx q[2];
rz(-0.060001686) q[2];
rz(1.9289198) q[3];
sx q[3];
rz(-2.4433177) q[3];
sx q[3];
rz(2.3771299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5948828) q[0];
sx q[0];
rz(-1.8130274) q[0];
sx q[0];
rz(-2.5643964) q[0];
rz(-0.82950854) q[1];
sx q[1];
rz(-2.0894158) q[1];
sx q[1];
rz(2.3037691) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5853607) q[0];
sx q[0];
rz(-0.62667003) q[0];
sx q[0];
rz(3.0100432) q[0];
rz(-pi) q[1];
rz(-0.98060645) q[2];
sx q[2];
rz(-1.6211121) q[2];
sx q[2];
rz(0.8129763) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7029757) q[1];
sx q[1];
rz(-1.9432406) q[1];
sx q[1];
rz(2.245642) q[1];
x q[2];
rz(-0.44562037) q[3];
sx q[3];
rz(-0.94404781) q[3];
sx q[3];
rz(-0.25139755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.025658) q[2];
sx q[2];
rz(-1.6178774) q[2];
sx q[2];
rz(-1.2238097) q[2];
rz(-2.6754248) q[3];
sx q[3];
rz(-2.7397082) q[3];
sx q[3];
rz(-0.60552067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8849628) q[0];
sx q[0];
rz(-1.4361359) q[0];
sx q[0];
rz(0.10109854) q[0];
rz(-0.413232) q[1];
sx q[1];
rz(-2.7006472) q[1];
sx q[1];
rz(-1.0879999) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4306385) q[0];
sx q[0];
rz(-1.0885518) q[0];
sx q[0];
rz(0.37094613) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3067683) q[2];
sx q[2];
rz(-1.2918279) q[2];
sx q[2];
rz(2.5503412) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.48641018) q[1];
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
rz(-pi) q[1];
rz(2.2130412) q[2];
sx q[2];
rz(-0.88625675) q[2];
sx q[2];
rz(-1.5555596) q[2];
rz(2.0689615) q[3];
sx q[3];
rz(-1.6173247) q[3];
sx q[3];
rz(-0.13256375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9685386) q[0];
sx q[0];
rz(-0.88468495) q[0];
sx q[0];
rz(-1.435085) q[0];
rz(-2.3834719) q[1];
sx q[1];
rz(-0.70223141) q[1];
sx q[1];
rz(0.10163669) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8827688) q[0];
sx q[0];
rz(-2.2166435) q[0];
sx q[0];
rz(0.96667883) q[0];
rz(-pi) q[1];
x q[1];
rz(0.36113895) q[2];
sx q[2];
rz(-1.1419347) q[2];
sx q[2];
rz(2.7782235) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1730301) q[1];
sx q[1];
rz(-0.30245879) q[1];
sx q[1];
rz(0.36424251) q[1];
rz(0.89438458) q[3];
sx q[3];
rz(-0.75856388) q[3];
sx q[3];
rz(-1.0762843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8018084) q[2];
sx q[2];
rz(-2.0399751) q[2];
sx q[2];
rz(1.8088809) q[2];
rz(-1.0821651) q[3];
sx q[3];
rz(-2.3887631) q[3];
sx q[3];
rz(1.4235206) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69671714) q[0];
sx q[0];
rz(-1.8901905) q[0];
sx q[0];
rz(-2.3984997) q[0];
rz(-0.7630868) q[1];
sx q[1];
rz(-0.63947314) q[1];
sx q[1];
rz(2.2230164) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42264274) q[0];
sx q[0];
rz(-3.115359) q[0];
sx q[0];
rz(-1.2573186) q[0];
x q[1];
rz(-0.95942504) q[2];
sx q[2];
rz(-0.55771962) q[2];
sx q[2];
rz(3.0447247) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2661675) q[1];
sx q[1];
rz(-1.2670749) q[1];
sx q[1];
rz(-2.6145302) q[1];
rz(-3.0094163) q[3];
sx q[3];
rz(-2.4130947) q[3];
sx q[3];
rz(-0.26154172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.7243728) q[2];
sx q[2];
rz(-1.959274) q[2];
sx q[2];
rz(2.7057538) q[2];
rz(-3.0271652) q[3];
sx q[3];
rz(-2.8947713) q[3];
sx q[3];
rz(2.54336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.8082387) q[0];
sx q[0];
rz(-0.55753189) q[0];
sx q[0];
rz(0.58407855) q[0];
rz(0.43560478) q[1];
sx q[1];
rz(-1.6807115) q[1];
sx q[1];
rz(-1.6216507) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7652033) q[0];
sx q[0];
rz(-2.0302773) q[0];
sx q[0];
rz(-2.8110519) q[0];
x q[1];
rz(-1.9797835) q[2];
sx q[2];
rz(-1.7715381) q[2];
sx q[2];
rz(-2.696352) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.73608172) q[1];
sx q[1];
rz(-2.9353432) q[1];
sx q[1];
rz(-3.095605) q[1];
rz(2.8826113) q[3];
sx q[3];
rz(-1.4121488) q[3];
sx q[3];
rz(0.64768956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.074177563) q[2];
sx q[2];
rz(-1.3520853) q[2];
sx q[2];
rz(2.0254859) q[2];
rz(1.3793147) q[3];
sx q[3];
rz(-1.1032871) q[3];
sx q[3];
rz(3.0232159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8922358) q[0];
sx q[0];
rz(-3.0638969) q[0];
sx q[0];
rz(-0.78999162) q[0];
rz(-3.0780011) q[1];
sx q[1];
rz(-0.60706943) q[1];
sx q[1];
rz(0.44073179) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9960956) q[0];
sx q[0];
rz(-0.4410559) q[0];
sx q[0];
rz(0.8121374) q[0];
rz(-pi) q[1];
rz(2.646801) q[2];
sx q[2];
rz(-1.2107009) q[2];
sx q[2];
rz(-0.39967135) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1524618) q[1];
sx q[1];
rz(-2.6351982) q[1];
sx q[1];
rz(-1.9146862) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3578254) q[3];
sx q[3];
rz(-2.9479288) q[3];
sx q[3];
rz(1.97399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.26238394) q[2];
sx q[2];
rz(-2.3449506) q[2];
sx q[2];
rz(2.4764496) q[2];
rz(2.3696259) q[3];
sx q[3];
rz(-0.77163458) q[3];
sx q[3];
rz(1.6316679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36122286) q[0];
sx q[0];
rz(-2.4936115) q[0];
sx q[0];
rz(-2.1369456) q[0];
rz(-1.4077582) q[1];
sx q[1];
rz(-1.6561457) q[1];
sx q[1];
rz(-1.8515324) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4125975) q[0];
sx q[0];
rz(-0.83461232) q[0];
sx q[0];
rz(-0.33552977) q[0];
rz(0.20310852) q[2];
sx q[2];
rz(-0.98047335) q[2];
sx q[2];
rz(-3.1121569) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2659246) q[1];
sx q[1];
rz(-0.34952119) q[1];
sx q[1];
rz(-0.65903475) q[1];
rz(1.8220002) q[3];
sx q[3];
rz(-0.61579865) q[3];
sx q[3];
rz(-1.3746266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1349692) q[2];
sx q[2];
rz(-1.1639736) q[2];
sx q[2];
rz(0.25035614) q[2];
rz(0.97744673) q[3];
sx q[3];
rz(-1.9340632) q[3];
sx q[3];
rz(1.4704963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1866495) q[0];
sx q[0];
rz(-1.432812) q[0];
sx q[0];
rz(-1.1695255) q[0];
rz(-2.3969338) q[1];
sx q[1];
rz(-0.79569334) q[1];
sx q[1];
rz(-2.4462499) q[1];
rz(1.5581239) q[2];
sx q[2];
rz(-1.4274393) q[2];
sx q[2];
rz(1.9096331) q[2];
rz(2.0674191) q[3];
sx q[3];
rz(-1.9575098) q[3];
sx q[3];
rz(-2.248275) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
