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
rz(-1.42000138759613) q[0];
sx q[0];
rz(5.32372585137422) q[0];
sx q[0];
rz(11.525679087631) q[0];
rz(-0.338484406471252) q[1];
sx q[1];
rz(1.65226331551606) q[1];
sx q[1];
rz(10.8545840740125) q[1];
cx q[1],q[0];
rz(-1.3022792339325) q[0];
sx q[0];
rz(2.84446180065209) q[0];
sx q[0];
rz(10.5797263145368) q[0];
rz(1.21665418148041) q[2];
sx q[2];
rz(2.14321890671784) q[2];
sx q[2];
rz(7.15711352824374) q[2];
cx q[2],q[1];
rz(-0.806244552135468) q[1];
sx q[1];
rz(4.99392155011232) q[1];
sx q[1];
rz(10.2608785390775) q[1];
rz(-2.87251758575439) q[3];
sx q[3];
rz(5.26139751275117) q[3];
sx q[3];
rz(10.8055635452191) q[3];
cx q[3],q[2];
rz(-1.13749837875366) q[2];
sx q[2];
rz(4.85529306729371) q[2];
sx q[2];
rz(12.6155404806058) q[2];
rz(2.5947425365448) q[3];
sx q[3];
rz(2.28737357457215) q[3];
sx q[3];
rz(8.12982604502841) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.262867242097855) q[0];
sx q[0];
rz(3.99223634799058) q[0];
sx q[0];
rz(9.38531502186462) q[0];
rz(0.706116855144501) q[1];
sx q[1];
rz(5.10895791848237) q[1];
sx q[1];
rz(11.3058472633283) q[1];
cx q[1],q[0];
rz(3.47760343551636) q[0];
sx q[0];
rz(5.12973562081391) q[0];
sx q[0];
rz(11.0593749046247) q[0];
rz(2.23725247383118) q[2];
sx q[2];
rz(5.15161577065522) q[2];
sx q[2];
rz(10.5365267753522) q[2];
cx q[2],q[1];
rz(-2.40521550178528) q[1];
sx q[1];
rz(4.74649730523164) q[1];
sx q[1];
rz(12.6085576772611) q[1];
rz(0.0665524080395699) q[3];
sx q[3];
rz(4.51887551148469) q[3];
sx q[3];
rz(8.96812704800769) q[3];
cx q[3],q[2];
rz(1.48271059989929) q[2];
sx q[2];
rz(5.27791729767854) q[2];
sx q[2];
rz(8.5073047041814) q[2];
rz(0.422580361366272) q[3];
sx q[3];
rz(5.03404334385926) q[3];
sx q[3];
rz(7.321170783035) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.680404007434845) q[0];
sx q[0];
rz(3.69643423159654) q[0];
sx q[0];
rz(10.3741606235425) q[0];
rz(1.89349114894867) q[1];
sx q[1];
rz(5.69948688347871) q[1];
sx q[1];
rz(7.79702053069278) q[1];
cx q[1],q[0];
rz(-0.682530522346497) q[0];
sx q[0];
rz(3.58624738653237) q[0];
sx q[0];
rz(11.594551062576) q[0];
rz(1.78691935539246) q[2];
sx q[2];
rz(5.18654933770234) q[2];
sx q[2];
rz(9.02468798159763) q[2];
cx q[2],q[1];
rz(-2.03841018676758) q[1];
sx q[1];
rz(5.00005451043183) q[1];
sx q[1];
rz(12.6436846017758) q[1];
rz(-0.413930028676987) q[3];
sx q[3];
rz(5.96998205979402) q[3];
sx q[3];
rz(11.8581712007444) q[3];
cx q[3],q[2];
rz(0.440243989229202) q[2];
sx q[2];
rz(2.01557186444337) q[2];
sx q[2];
rz(11.0429421424787) q[2];
rz(2.82738876342773) q[3];
sx q[3];
rz(5.25755563576753) q[3];
sx q[3];
rz(7.92269775866672) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.03581285476685) q[0];
sx q[0];
rz(3.33296027977998) q[0];
sx q[0];
rz(9.58187035321399) q[0];
rz(3.00639033317566) q[1];
sx q[1];
rz(3.92670026619966) q[1];
sx q[1];
rz(9.0915579855363) q[1];
cx q[1],q[0];
rz(-1.55578088760376) q[0];
sx q[0];
rz(5.58526483376557) q[0];
sx q[0];
rz(11.7410716772) q[0];
rz(-2.9178102016449) q[2];
sx q[2];
rz(5.54229679902131) q[2];
sx q[2];
rz(6.51349470614597) q[2];
cx q[2],q[1];
rz(3.13184785842896) q[1];
sx q[1];
rz(1.83418432076509) q[1];
sx q[1];
rz(3.5034241437833) q[1];
rz(-0.35857692360878) q[3];
sx q[3];
rz(5.53809824784333) q[3];
sx q[3];
rz(11.6803750753324) q[3];
cx q[3],q[2];
rz(0.80804431438446) q[2];
sx q[2];
rz(2.49258884985978) q[2];
sx q[2];
rz(7.9760700225751) q[2];
rz(2.40888047218323) q[3];
sx q[3];
rz(5.06356659730012) q[3];
sx q[3];
rz(4.75134370326206) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-2.69199585914612) q[0];
sx q[0];
rz(4.62822547753388) q[0];
sx q[0];
rz(9.65633348225757) q[0];
rz(0.721365034580231) q[1];
sx q[1];
rz(4.03272745211656) q[1];
sx q[1];
rz(10.2732073426168) q[1];
cx q[1],q[0];
rz(1.41302943229675) q[0];
sx q[0];
rz(3.74429753621156) q[0];
sx q[0];
rz(8.80828932522937) q[0];
rz(0.404506385326385) q[2];
sx q[2];
rz(5.21878138382966) q[2];
sx q[2];
rz(9.60410145520374) q[2];
cx q[2],q[1];
rz(0.906815469264984) q[1];
sx q[1];
rz(3.53580010135705) q[1];
sx q[1];
rz(11.1787062644879) q[1];
rz(-4.30013608932495) q[3];
sx q[3];
rz(0.257671507196971) q[3];
sx q[3];
rz(10.1098750591199) q[3];
cx q[3],q[2];
rz(-1.1344176530838) q[2];
sx q[2];
rz(3.7448132952028) q[2];
sx q[2];
rz(11.3654320001523) q[2];
rz(0.689315557479858) q[3];
sx q[3];
rz(5.11316195328767) q[3];
sx q[3];
rz(9.9905040025632) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.98179626464844) q[0];
sx q[0];
rz(1.60982755024964) q[0];
sx q[0];
rz(11.2181665658872) q[0];
rz(-2.42344045639038) q[1];
sx q[1];
rz(3.71287957032258) q[1];
sx q[1];
rz(11.274615263931) q[1];
cx q[1],q[0];
rz(-0.517070889472961) q[0];
sx q[0];
rz(5.87070575554902) q[0];
sx q[0];
rz(8.44256517886325) q[0];
rz(3.74019479751587) q[2];
sx q[2];
rz(4.39843240578706) q[2];
sx q[2];
rz(10.4281838893811) q[2];
cx q[2],q[1];
rz(0.764522194862366) q[1];
sx q[1];
rz(1.65652802784974) q[1];
sx q[1];
rz(8.02315638064548) q[1];
rz(1.51025283336639) q[3];
sx q[3];
rz(5.88715663750703) q[3];
sx q[3];
rz(10.1939216017644) q[3];
cx q[3],q[2];
rz(1.98880875110626) q[2];
sx q[2];
rz(5.31816426117951) q[2];
sx q[2];
rz(9.82688609360858) q[2];
rz(2.13573384284973) q[3];
sx q[3];
rz(3.36496027012403) q[3];
sx q[3];
rz(7.98739538191959) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.417563676834106) q[0];
sx q[0];
rz(3.28983131249482) q[0];
sx q[0];
rz(10.754368042938) q[0];
rz(1.81364643573761) q[1];
sx q[1];
rz(1.41683939297731) q[1];
sx q[1];
rz(6.73483107089206) q[1];
cx q[1],q[0];
rz(0.992579162120819) q[0];
sx q[0];
rz(4.27044430573518) q[0];
sx q[0];
rz(9.32391079365417) q[0];
rz(1.99597859382629) q[2];
sx q[2];
rz(1.69847920735414) q[2];
sx q[2];
rz(9.04421541690036) q[2];
cx q[2],q[1];
rz(1.82800805568695) q[1];
sx q[1];
rz(0.196889551477977) q[1];
sx q[1];
rz(6.02189157008334) q[1];
rz(-2.3673243522644) q[3];
sx q[3];
rz(2.30832705100114) q[3];
sx q[3];
rz(9.31302261947795) q[3];
cx q[3],q[2];
rz(0.202465012669563) q[2];
sx q[2];
rz(4.51345280011232) q[2];
sx q[2];
rz(9.71264923214122) q[2];
rz(1.48160827159882) q[3];
sx q[3];
rz(3.99325993855531) q[3];
sx q[3];
rz(10.8816359996717) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.579563677310944) q[0];
sx q[0];
rz(4.10211143096025) q[0];
sx q[0];
rz(8.87185958623096) q[0];
rz(-1.25060725212097) q[1];
sx q[1];
rz(3.55184796650941) q[1];
sx q[1];
rz(8.04432067870303) q[1];
cx q[1],q[0];
rz(-0.569262146949768) q[0];
sx q[0];
rz(2.49118605454499) q[0];
sx q[0];
rz(7.22800586222812) q[0];
rz(2.55437183380127) q[2];
sx q[2];
rz(4.81382778485353) q[2];
sx q[2];
rz(7.59562501906558) q[2];
cx q[2],q[1];
rz(-0.985017240047455) q[1];
sx q[1];
rz(3.38216552336747) q[1];
sx q[1];
rz(13.9634461164395) q[1];
rz(0.172354802489281) q[3];
sx q[3];
rz(5.50561753113801) q[3];
sx q[3];
rz(8.94320300816699) q[3];
cx q[3],q[2];
rz(0.550142526626587) q[2];
sx q[2];
rz(4.93530610402162) q[2];
sx q[2];
rz(6.99828622340366) q[2];
rz(-0.0446557849645615) q[3];
sx q[3];
rz(2.55845925410325) q[3];
sx q[3];
rz(11.9379899263303) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.1735246181488) q[0];
sx q[0];
rz(5.5139423926645) q[0];
sx q[0];
rz(8.78389099835559) q[0];
rz(1.64546537399292) q[1];
sx q[1];
rz(3.5265453179651) q[1];
sx q[1];
rz(8.77326039075061) q[1];
cx q[1],q[0];
rz(-4.34234952926636) q[0];
sx q[0];
rz(4.02389207680757) q[0];
sx q[0];
rz(10.6606454610746) q[0];
rz(-0.233853057026863) q[2];
sx q[2];
rz(7.37557998498017) q[2];
sx q[2];
rz(12.6209485292356) q[2];
cx q[2],q[1];
rz(5.289466381073) q[1];
sx q[1];
rz(5.35618820984895) q[1];
sx q[1];
rz(10.5293291568677) q[1];
rz(0.348179489374161) q[3];
sx q[3];
rz(4.89236012299592) q[3];
sx q[3];
rz(8.9090238571088) q[3];
cx q[3],q[2];
rz(-1.73518252372742) q[2];
sx q[2];
rz(2.09690216382081) q[2];
sx q[2];
rz(12.1082885026853) q[2];
rz(-0.0118907513096929) q[3];
sx q[3];
rz(4.20227483113343) q[3];
sx q[3];
rz(9.40301214008733) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.402820140123367) q[0];
sx q[0];
rz(3.33644765813882) q[0];
sx q[0];
rz(9.4537549059759) q[0];
rz(-0.926523566246033) q[1];
sx q[1];
rz(4.59007528622682) q[1];
sx q[1];
rz(9.01851857303783) q[1];
cx q[1],q[0];
rz(0.479325205087662) q[0];
sx q[0];
rz(2.86135840614373) q[0];
sx q[0];
rz(11.026666378967) q[0];
rz(0.455300718545914) q[2];
sx q[2];
rz(5.58497801621492) q[2];
sx q[2];
rz(10.355744457237) q[2];
cx q[2],q[1];
rz(-4.09029626846313) q[1];
sx q[1];
rz(1.27928415139253) q[1];
sx q[1];
rz(7.88647637366458) q[1];
rz(-0.77099734544754) q[3];
sx q[3];
rz(2.16350284417207) q[3];
sx q[3];
rz(10.5143804311673) q[3];
cx q[3],q[2];
rz(2.37090921401978) q[2];
sx q[2];
rz(4.63969722588594) q[2];
sx q[2];
rz(11.9715978860776) q[2];
rz(2.4316782951355) q[3];
sx q[3];
rz(4.09865215619142) q[3];
sx q[3];
rz(8.47944215535327) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.58890771865845) q[0];
sx q[0];
rz(5.38059178193147) q[0];
sx q[0];
rz(10.2519810557286) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(1.20811879634857) q[1];
sx q[1];
rz(1.79580095608766) q[1];
sx q[1];
rz(9.67617369293376) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(-0.571115374565125) q[2];
sx q[2];
rz(4.21943715413148) q[2];
sx q[2];
rz(9.87710968255206) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(-0.609676539897919) q[3];
sx q[3];
rz(3.57537654240663) q[3];
sx q[3];
rz(11.3479261159818) q[3];
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
