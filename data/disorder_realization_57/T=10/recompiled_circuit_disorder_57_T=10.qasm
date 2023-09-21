OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg meas[4];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.527430713176727) q[0];
sx q[0];
rz(3.95137402613694) q[0];
sx q[0];
rz(9.95617190598651) q[0];
rz(0.237513795495033) q[1];
sx q[1];
rz(1.3637434562021) q[1];
sx q[1];
rz(10.6278282165448) q[1];
cx q[1],q[0];
rz(4.99942398071289) q[0];
sx q[0];
rz(5.99537483056123) q[0];
sx q[0];
rz(10.1596927404325) q[0];
rz(-0.216993987560272) q[2];
sx q[2];
rz(1.81976822217042) q[2];
sx q[2];
rz(8.70017549990817) q[2];
cx q[2],q[1];
rz(-0.141292780637741) q[1];
sx q[1];
rz(4.94718268712098) q[1];
sx q[1];
rz(12.5077862501065) q[1];
rz(-1.86373782157898) q[3];
sx q[3];
rz(4.74786189396913) q[3];
sx q[3];
rz(12.002629017822) q[3];
cx q[3],q[2];
rz(0.276658594608307) q[2];
sx q[2];
rz(4.27327159245545) q[2];
sx q[2];
rz(9.62386185526058) q[2];
rz(1.52786004543304) q[3];
sx q[3];
rz(3.63854232628877) q[3];
sx q[3];
rz(8.69069603680774) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.14282250404358) q[0];
sx q[0];
rz(3.26876636047895) q[0];
sx q[0];
rz(8.60536805390521) q[0];
rz(-0.285851299762726) q[1];
sx q[1];
rz(4.05552301009233) q[1];
sx q[1];
rz(11.2999068260114) q[1];
cx q[1],q[0];
rz(-1.09853839874268) q[0];
sx q[0];
rz(2.17159542639787) q[0];
sx q[0];
rz(13.0310544729154) q[0];
rz(-0.404475003480911) q[2];
sx q[2];
rz(1.35044875939424) q[2];
sx q[2];
rz(7.88274130820438) q[2];
cx q[2],q[1];
rz(3.75209879875183) q[1];
sx q[1];
rz(6.83581987221772) q[1];
sx q[1];
rz(11.8289930581967) q[1];
rz(0.451617509126663) q[3];
sx q[3];
rz(1.04188099701936) q[3];
sx q[3];
rz(12.43922326564) q[3];
cx q[3],q[2];
rz(1.15660285949707) q[2];
sx q[2];
rz(5.06982508500154) q[2];
sx q[2];
rz(7.445154285423) q[2];
rz(3.05442881584167) q[3];
sx q[3];
rz(1.63790837128694) q[3];
sx q[3];
rz(9.4986869379799) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.53482985496521) q[0];
sx q[0];
rz(4.26229527791078) q[0];
sx q[0];
rz(9.12117750047847) q[0];
rz(-1.7595694065094) q[1];
sx q[1];
rz(4.41777387459809) q[1];
sx q[1];
rz(13.2138707399289) q[1];
cx q[1],q[0];
rz(-0.305953681468964) q[0];
sx q[0];
rz(1.52579716046388) q[0];
sx q[0];
rz(10.3609844803731) q[0];
rz(-0.961704015731812) q[2];
sx q[2];
rz(3.95590809186036) q[2];
sx q[2];
rz(8.98440582155391) q[2];
cx q[2],q[1];
rz(1.46090471744537) q[1];
sx q[1];
rz(3.74244228203828) q[1];
sx q[1];
rz(11.8851096391599) q[1];
rz(0.697528719902039) q[3];
sx q[3];
rz(4.23444917996461) q[3];
sx q[3];
rz(11.5124895334165) q[3];
cx q[3],q[2];
rz(1.29937624931335) q[2];
sx q[2];
rz(2.35895779927308) q[2];
sx q[2];
rz(10.9856388330381) q[2];
rz(2.1598699092865) q[3];
sx q[3];
rz(5.00365463097627) q[3];
sx q[3];
rz(10.0681899547498) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.90755391120911) q[0];
sx q[0];
rz(2.30383947690065) q[0];
sx q[0];
rz(10.2966762542646) q[0];
rz(-2.75432276725769) q[1];
sx q[1];
rz(2.49881199200685) q[1];
sx q[1];
rz(12.2753247976224) q[1];
cx q[1],q[0];
rz(-0.298453330993652) q[0];
sx q[0];
rz(1.79945865471894) q[0];
sx q[0];
rz(8.83295325040027) q[0];
rz(-1.86371517181396) q[2];
sx q[2];
rz(4.40410903294618) q[2];
sx q[2];
rz(7.04702327250644) q[2];
cx q[2],q[1];
rz(3.31268310546875) q[1];
sx q[1];
rz(5.39576783974702) q[1];
sx q[1];
rz(10.3239345908086) q[1];
rz(0.244922280311584) q[3];
sx q[3];
rz(3.73186001380021) q[3];
sx q[3];
rz(9.31065416931316) q[3];
cx q[3],q[2];
rz(1.77240931987762) q[2];
sx q[2];
rz(5.48824778397615) q[2];
sx q[2];
rz(8.88072142600223) q[2];
rz(0.509285151958466) q[3];
sx q[3];
rz(5.21930948098237) q[3];
sx q[3];
rz(10.2907519698064) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.139570251107216) q[0];
sx q[0];
rz(2.11525705655152) q[0];
sx q[0];
rz(10.0698445796888) q[0];
rz(0.625726103782654) q[1];
sx q[1];
rz(4.36521700223024) q[1];
sx q[1];
rz(10.0549567103307) q[1];
cx q[1],q[0];
rz(1.52983129024506) q[0];
sx q[0];
rz(5.28695455391938) q[0];
sx q[0];
rz(9.2792237162511) q[0];
rz(3.89242196083069) q[2];
sx q[2];
rz(4.41725650628144) q[2];
sx q[2];
rz(5.59460661410495) q[2];
cx q[2],q[1];
rz(0.678179204463959) q[1];
sx q[1];
rz(5.6971801837259) q[1];
sx q[1];
rz(9.83739343880817) q[1];
rz(1.8120368719101) q[3];
sx q[3];
rz(2.57062283356721) q[3];
sx q[3];
rz(6.0007374048154) q[3];
cx q[3],q[2];
rz(-0.753829181194305) q[2];
sx q[2];
rz(1.32396534283692) q[2];
sx q[2];
rz(10.7889443397443) q[2];
rz(-1.24986135959625) q[3];
sx q[3];
rz(4.35721638997132) q[3];
sx q[3];
rz(8.50331923960849) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.08677780628204) q[0];
sx q[0];
rz(3.7253037412935) q[0];
sx q[0];
rz(8.70794997214481) q[0];
rz(-3.57930541038513) q[1];
sx q[1];
rz(3.82203939755494) q[1];
sx q[1];
rz(11.752669787399) q[1];
cx q[1],q[0];
rz(-0.871562242507935) q[0];
sx q[0];
rz(3.36420747836167) q[0];
sx q[0];
rz(9.71280739306613) q[0];
rz(0.93136602640152) q[2];
sx q[2];
rz(5.45139470894868) q[2];
sx q[2];
rz(9.12162808179065) q[2];
cx q[2],q[1];
rz(-2.59722828865051) q[1];
sx q[1];
rz(4.37889102299745) q[1];
sx q[1];
rz(11.7935564279477) q[1];
rz(-0.709909856319427) q[3];
sx q[3];
rz(6.73618808587129) q[3];
sx q[3];
rz(8.59121069907352) q[3];
cx q[3],q[2];
rz(1.25632357597351) q[2];
sx q[2];
rz(3.57907390792901) q[2];
sx q[2];
rz(8.23711857794925) q[2];
rz(-2.20962834358215) q[3];
sx q[3];
rz(4.27560511429841) q[3];
sx q[3];
rz(9.79595342873737) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.15357828140259) q[0];
sx q[0];
rz(4.17052117188508) q[0];
sx q[0];
rz(10.1108614563863) q[0];
rz(1.52308058738708) q[1];
sx q[1];
rz(4.11580357153947) q[1];
sx q[1];
rz(8.76151744126483) q[1];
cx q[1],q[0];
rz(-1.23847568035126) q[0];
sx q[0];
rz(3.79229328234727) q[0];
sx q[0];
rz(9.38530262782379) q[0];
rz(1.86134672164917) q[2];
sx q[2];
rz(3.79972639878327) q[2];
sx q[2];
rz(10.9550416231076) q[2];
cx q[2],q[1];
rz(2.75463318824768) q[1];
sx q[1];
rz(2.61357507308061) q[1];
sx q[1];
rz(9.01628593205615) q[1];
rz(2.25606298446655) q[3];
sx q[3];
rz(3.82420191367204) q[3];
sx q[3];
rz(9.63519174455806) q[3];
cx q[3],q[2];
rz(2.47458100318909) q[2];
sx q[2];
rz(3.97908535798127) q[2];
sx q[2];
rz(9.43985826186045) q[2];
rz(1.01174259185791) q[3];
sx q[3];
rz(5.16705432732637) q[3];
sx q[3];
rz(9.99619261025592) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.0529772341251373) q[0];
sx q[0];
rz(0.732607277231761) q[0];
sx q[0];
rz(9.3173892557542) q[0];
rz(-0.304741621017456) q[1];
sx q[1];
rz(4.46008232434327) q[1];
sx q[1];
rz(9.8825565636079) q[1];
cx q[1],q[0];
rz(0.251541972160339) q[0];
sx q[0];
rz(5.5596708377176) q[0];
sx q[0];
rz(9.92199382781192) q[0];
rz(-1.00502145290375) q[2];
sx q[2];
rz(3.97275576193864) q[2];
sx q[2];
rz(10.0891454577367) q[2];
cx q[2],q[1];
rz(-0.0195518601685762) q[1];
sx q[1];
rz(4.47960141499574) q[1];
sx q[1];
rz(12.2003135442655) q[1];
rz(1.02022445201874) q[3];
sx q[3];
rz(2.51591608126695) q[3];
sx q[3];
rz(10.953233575813) q[3];
cx q[3],q[2];
rz(-2.77659392356873) q[2];
sx q[2];
rz(4.59805432160432) q[2];
sx q[2];
rz(7.71467051505252) q[2];
rz(1.71639037132263) q[3];
sx q[3];
rz(3.76834985812242) q[3];
sx q[3];
rz(8.13858506678745) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.154370561242104) q[0];
sx q[0];
rz(2.66193375189836) q[0];
sx q[0];
rz(11.6128837823789) q[0];
rz(1.32572388648987) q[1];
sx q[1];
rz(3.63453519542749) q[1];
sx q[1];
rz(10.0095172882001) q[1];
cx q[1],q[0];
rz(0.0277715865522623) q[0];
sx q[0];
rz(2.63915327389772) q[0];
sx q[0];
rz(10.7611760854642) q[0];
rz(3.50727128982544) q[2];
sx q[2];
rz(1.0569910128885) q[2];
sx q[2];
rz(9.9300324678342) q[2];
cx q[2],q[1];
rz(0.194962754845619) q[1];
sx q[1];
rz(1.13190356095368) q[1];
sx q[1];
rz(9.50131413935825) q[1];
rz(0.483804881572723) q[3];
sx q[3];
rz(3.44171059330041) q[3];
sx q[3];
rz(9.12394646405383) q[3];
cx q[3],q[2];
rz(1.33283019065857) q[2];
sx q[2];
rz(2.32138827641542) q[2];
sx q[2];
rz(9.70104832052394) q[2];
rz(0.551094651222229) q[3];
sx q[3];
rz(4.88376704056794) q[3];
sx q[3];
rz(9.04111102818652) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.06934142112732) q[0];
sx q[0];
rz(3.65579363902146) q[0];
sx q[0];
rz(10.0802615046422) q[0];
rz(1.15707015991211) q[1];
sx q[1];
rz(4.90166309674317) q[1];
sx q[1];
rz(11.0438372850339) q[1];
cx q[1],q[0];
rz(-0.0392783954739571) q[0];
sx q[0];
rz(3.88003930647904) q[0];
sx q[0];
rz(10.5260152578275) q[0];
rz(2.25879144668579) q[2];
sx q[2];
rz(1.30324593384797) q[2];
sx q[2];
rz(6.88163850306674) q[2];
cx q[2],q[1];
rz(1.57877051830292) q[1];
sx q[1];
rz(2.49145034153993) q[1];
sx q[1];
rz(11.269195175163) q[1];
rz(-0.811509847640991) q[3];
sx q[3];
rz(5.72890296776826) q[3];
sx q[3];
rz(10.049393272392) q[3];
cx q[3],q[2];
rz(2.98605012893677) q[2];
sx q[2];
rz(3.54938554962213) q[2];
sx q[2];
rz(10.2669196486394) q[2];
rz(1.53676736354828) q[3];
sx q[3];
rz(1.38048806984956) q[3];
sx q[3];
rz(9.39824903420314) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.29795169830322) q[0];
sx q[0];
rz(1.15402963955934) q[0];
sx q[0];
rz(9.00463700889751) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(-0.783214449882507) q[1];
sx q[1];
rz(3.83677438099916) q[1];
sx q[1];
rz(11.0421257972638) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(1.51700675487518) q[2];
sx q[2];
rz(3.35136556823785) q[2];
sx q[2];
rz(14.2158360242765) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(-0.363568156957626) q[3];
sx q[3];
rz(6.21842852433259) q[3];
sx q[3];
rz(8.18063197135135) q[3];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[2] -> meas[2];
measure q[3] -> meas[3];