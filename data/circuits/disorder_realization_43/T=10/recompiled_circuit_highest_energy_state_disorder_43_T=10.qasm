OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg meas[4];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.690052032470703) q[0];
sx q[0];
rz(8.57418409188325) q[0];
sx q[0];
rz(9.22467145919009) q[0];
rz(-0.198356807231903) q[1];
sx q[1];
rz(5.76016465027864) q[1];
sx q[1];
rz(10.7566071510236) q[1];
cx q[1],q[0];
rz(-5.09766674041748) q[0];
sx q[0];
rz(8.08543458779389) q[0];
sx q[0];
rz(8.9737723827283) q[0];
rz(3.63177490234375) q[2];
sx q[2];
rz(1.02875438530976) q[2];
sx q[2];
rz(4.4619669675748) q[2];
cx q[2],q[1];
rz(-1.68649041652679) q[1];
sx q[1];
rz(5.38977924187715) q[1];
sx q[1];
rz(12.6280748605649) q[1];
rz(0.159456342458725) q[3];
sx q[3];
rz(6.7644089778238) q[3];
sx q[3];
rz(6.38336822985812) q[3];
cx q[3],q[2];
rz(3.06938219070435) q[2];
sx q[2];
rz(3.50241530139978) q[2];
sx q[2];
rz(4.45617578028842) q[2];
rz(-2.23831105232239) q[3];
sx q[3];
rz(1.77777806122834) q[3];
sx q[3];
rz(12.2204632520597) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.271799296140671) q[0];
sx q[0];
rz(1.48577395279939) q[0];
sx q[0];
rz(10.2805697083394) q[0];
rz(-0.774041473865509) q[1];
sx q[1];
rz(4.30553320248658) q[1];
sx q[1];
rz(11.7785076856534) q[1];
cx q[1],q[0];
rz(-2.6007399559021) q[0];
sx q[0];
rz(2.88982805808122) q[0];
sx q[0];
rz(6.19397375582858) q[0];
rz(3.18399381637573) q[2];
sx q[2];
rz(5.05633607705171) q[2];
sx q[2];
rz(16.2016787290494) q[2];
cx q[2],q[1];
rz(0.987092435359955) q[1];
sx q[1];
rz(-0.561141578359059) q[1];
sx q[1];
rz(8.74946866034671) q[1];
rz(0.90372371673584) q[3];
sx q[3];
rz(4.87822678883607) q[3];
sx q[3];
rz(7.47820959090396) q[3];
cx q[3],q[2];
rz(3.07561254501343) q[2];
sx q[2];
rz(8.1947113593393) q[2];
sx q[2];
rz(6.67928502558872) q[2];
rz(-1.57054495811462) q[3];
sx q[3];
rz(-1.50382646719878) q[3];
sx q[3];
rz(10.7802475452344) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.95392549037933) q[0];
sx q[0];
rz(4.7879697402292) q[0];
sx q[0];
rz(9.48538806884691) q[0];
rz(2.45687794685364) q[1];
sx q[1];
rz(4.54185334046418) q[1];
sx q[1];
rz(12.2273247003476) q[1];
cx q[1],q[0];
rz(-3.54450964927673) q[0];
sx q[0];
rz(3.72501304944093) q[0];
sx q[0];
rz(7.81636891364261) q[0];
rz(8.42842960357666) q[2];
sx q[2];
rz(2.91746047337586) q[2];
sx q[2];
rz(14.0382771253507) q[2];
cx q[2],q[1];
rz(-0.782796919345856) q[1];
sx q[1];
rz(2.72272449930245) q[1];
sx q[1];
rz(11.8435163259427) q[1];
rz(0.0117118917405605) q[3];
sx q[3];
rz(5.48918929894502) q[3];
sx q[3];
rz(6.65028140544101) q[3];
cx q[3],q[2];
rz(2.65913486480713) q[2];
sx q[2];
rz(1.86379161675508) q[2];
sx q[2];
rz(13.0368320703427) q[2];
rz(-0.862369239330292) q[3];
sx q[3];
rz(2.69263977010781) q[3];
sx q[3];
rz(7.86320028304263) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.77931332588196) q[0];
sx q[0];
rz(5.16474214394624) q[0];
sx q[0];
rz(5.87516663073703) q[0];
rz(1.79040026664734) q[1];
sx q[1];
rz(4.98570957978303) q[1];
sx q[1];
rz(4.96195027827426) q[1];
cx q[1],q[0];
rz(3.35622549057007) q[0];
sx q[0];
rz(4.21186617215211) q[0];
sx q[0];
rz(15.3619794607083) q[0];
rz(-3.25339818000793) q[2];
sx q[2];
rz(6.94139877160127) q[2];
sx q[2];
rz(7.64615545272037) q[2];
cx q[2],q[1];
rz(6.01279783248901) q[1];
sx q[1];
rz(3.28318146069581) q[1];
sx q[1];
rz(7.3130616903226) q[1];
rz(1.14369904994965) q[3];
sx q[3];
rz(5.09988120396669) q[3];
sx q[3];
rz(10.3856413125913) q[3];
cx q[3],q[2];
rz(-4.36461019515991) q[2];
sx q[2];
rz(1.43179801304872) q[2];
sx q[2];
rz(11.8435165643613) q[2];
rz(0.849607884883881) q[3];
sx q[3];
rz(4.31416419346864) q[3];
sx q[3];
rz(10.4626487255017) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.95491135120392) q[0];
sx q[0];
rz(3.6584930737787) q[0];
sx q[0];
rz(6.84315202235385) q[0];
rz(-3.34666705131531) q[1];
sx q[1];
rz(5.0264951308542) q[1];
sx q[1];
rz(6.57427403926059) q[1];
cx q[1],q[0];
rz(-1.73228919506073) q[0];
sx q[0];
rz(1.83900526364381) q[0];
sx q[0];
rz(6.10755631922885) q[0];
rz(-1.5093549489975) q[2];
sx q[2];
rz(3.94179353316362) q[2];
sx q[2];
rz(9.07348615526363) q[2];
cx q[2],q[1];
rz(-6.41097211837769) q[1];
sx q[1];
rz(4.9108716567331) q[1];
sx q[1];
rz(16.4119562864225) q[1];
rz(0.160334825515747) q[3];
sx q[3];
rz(5.45780530770356) q[3];
sx q[3];
rz(10.8061188220899) q[3];
cx q[3],q[2];
rz(1.22289550304413) q[2];
sx q[2];
rz(4.61997845967347) q[2];
sx q[2];
rz(5.98488662242099) q[2];
rz(1.97226238250732) q[3];
sx q[3];
rz(4.10852125485475) q[3];
sx q[3];
rz(7.91351613997623) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.83396124839783) q[0];
sx q[0];
rz(4.60222485859925) q[0];
sx q[0];
rz(8.82237676381274) q[0];
rz(6.6547327041626) q[1];
sx q[1];
rz(4.68640783627565) q[1];
sx q[1];
rz(4.1733550786893) q[1];
cx q[1],q[0];
rz(2.27157974243164) q[0];
sx q[0];
rz(4.53813412983949) q[0];
sx q[0];
rz(10.9794952630918) q[0];
rz(-2.93390536308289) q[2];
sx q[2];
rz(4.31549206574494) q[2];
sx q[2];
rz(8.9924274444501) q[2];
cx q[2],q[1];
rz(0.183730155229568) q[1];
sx q[1];
rz(2.0477843602472) q[1];
sx q[1];
rz(12.3535666227262) q[1];
rz(0.647531151771545) q[3];
sx q[3];
rz(4.23491624196107) q[3];
sx q[3];
rz(11.060551261894) q[3];
cx q[3],q[2];
rz(-0.804566979408264) q[2];
sx q[2];
rz(4.24293104012544) q[2];
sx q[2];
rz(14.1892719030301) q[2];
rz(0.848414599895477) q[3];
sx q[3];
rz(4.01685932477052) q[3];
sx q[3];
rz(12.6059810876767) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.11617159843445) q[0];
sx q[0];
rz(3.09053675283725) q[0];
sx q[0];
rz(11.5684711694638) q[0];
rz(0.274670302867889) q[1];
sx q[1];
rz(4.02172860701615) q[1];
sx q[1];
rz(10.7956478357236) q[1];
cx q[1],q[0];
rz(1.27517282962799) q[0];
sx q[0];
rz(2.20179704030091) q[0];
sx q[0];
rz(12.5830271005551) q[0];
rz(1.72452318668365) q[2];
sx q[2];
rz(5.06444791157777) q[2];
sx q[2];
rz(11.0219039678495) q[2];
cx q[2],q[1];
rz(0.26105597615242) q[1];
sx q[1];
rz(2.47891006072099) q[1];
sx q[1];
rz(4.31463572978183) q[1];
rz(-0.771057188510895) q[3];
sx q[3];
rz(4.13660666544969) q[3];
sx q[3];
rz(6.76597783564731) q[3];
cx q[3],q[2];
rz(2.28445982933044) q[2];
sx q[2];
rz(0.942281158762523) q[2];
sx q[2];
rz(6.31283161639377) q[2];
rz(-0.615214109420776) q[3];
sx q[3];
rz(4.88636508782441) q[3];
sx q[3];
rz(9.76711595653697) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.0895594730973244) q[0];
sx q[0];
rz(5.82089868386323) q[0];
sx q[0];
rz(8.64172205924197) q[0];
rz(-1.1198593378067) q[1];
sx q[1];
rz(3.80870661337907) q[1];
sx q[1];
rz(10.8227523326795) q[1];
cx q[1],q[0];
rz(-0.0722948387265205) q[0];
sx q[0];
rz(5.1552404483133) q[0];
sx q[0];
rz(11.4504399061124) q[0];
rz(0.284420222043991) q[2];
sx q[2];
rz(4.4878900368982) q[2];
sx q[2];
rz(9.03628442286655) q[2];
cx q[2],q[1];
rz(0.840613424777985) q[1];
sx q[1];
rz(5.61524215539033) q[1];
sx q[1];
rz(9.43901852927312) q[1];
rz(-1.78461515903473) q[3];
sx q[3];
rz(5.48733869393403) q[3];
sx q[3];
rz(8.72382853030368) q[3];
cx q[3],q[2];
rz(0.295842587947845) q[2];
sx q[2];
rz(2.55518987973268) q[2];
sx q[2];
rz(10.7455156803052) q[2];
rz(-1.01975059509277) q[3];
sx q[3];
rz(4.56557885010774) q[3];
sx q[3];
rz(10.8446835040967) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.94186973571777) q[0];
sx q[0];
rz(2.04827347596223) q[0];
sx q[0];
rz(11.6119687318723) q[0];
rz(1.15425395965576) q[1];
sx q[1];
rz(5.63710919220979) q[1];
sx q[1];
rz(7.3676066160123) q[1];
cx q[1],q[0];
rz(-1.08237755298615) q[0];
sx q[0];
rz(3.85126927693421) q[0];
sx q[0];
rz(9.48952715694114) q[0];
rz(-0.280867516994476) q[2];
sx q[2];
rz(4.46675148804719) q[2];
sx q[2];
rz(3.5421786069791) q[2];
cx q[2],q[1];
rz(6.32818651199341) q[1];
sx q[1];
rz(3.61970743735368) q[1];
sx q[1];
rz(6.44085047244235) q[1];
rz(1.37280416488647) q[3];
sx q[3];
rz(5.52471295197541) q[3];
sx q[3];
rz(14.7713684797208) q[3];
cx q[3],q[2];
rz(0.0831009745597839) q[2];
sx q[2];
rz(6.47506085236604) q[2];
sx q[2];
rz(9.86189678906604) q[2];
rz(7.62656021118164) q[3];
sx q[3];
rz(4.49298134644563) q[3];
sx q[3];
rz(9.37940645068094) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-3.09923911094666) q[0];
sx q[0];
rz(4.14206114609773) q[0];
sx q[0];
rz(7.69371089934512) q[0];
rz(-0.220774158835411) q[1];
sx q[1];
rz(5.75194803078706) q[1];
sx q[1];
rz(8.49154726266071) q[1];
cx q[1],q[0];
rz(-1.27569031715393) q[0];
sx q[0];
rz(2.31610158284242) q[0];
sx q[0];
rz(9.44335288963422) q[0];
rz(3.50191926956177) q[2];
sx q[2];
rz(5.30012789567048) q[2];
sx q[2];
rz(9.08247954248592) q[2];
cx q[2],q[1];
rz(0.905518412590027) q[1];
sx q[1];
rz(2.33744403918321) q[1];
sx q[1];
rz(11.1797033309857) q[1];
rz(1.68269622325897) q[3];
sx q[3];
rz(0.941169412928172) q[3];
sx q[3];
rz(5.31795117854282) q[3];
cx q[3],q[2];
rz(3.85179734230042) q[2];
sx q[2];
rz(3.48381823499734) q[2];
sx q[2];
rz(6.9146504163663) q[2];
rz(-3.86059784889221) q[3];
sx q[3];
rz(1.34716657002503) q[3];
sx q[3];
rz(9.6662594884555) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.635425686836243) q[0];
sx q[0];
rz(0.271871956186839) q[0];
sx q[0];
rz(6.59062192439243) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(-3.80481815338135) q[1];
sx q[1];
rz(2.66223380168016) q[1];
sx q[1];
rz(15.1923918485562) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(3.41660523414612) q[2];
sx q[2];
rz(4.1703132708841) q[2];
sx q[2];
rz(5.83468625544711) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(1.94619464874268) q[3];
sx q[3];
rz(1.97070983250672) q[3];
sx q[3];
rz(15.0851902723233) q[3];
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
