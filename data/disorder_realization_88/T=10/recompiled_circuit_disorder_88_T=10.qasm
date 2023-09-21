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
rz(-0.247034385800362) q[0];
sx q[0];
rz(4.39727536042268) q[0];
sx q[0];
rz(9.75274047850772) q[0];
rz(-0.234499603509903) q[1];
sx q[1];
rz(3.34267153044278) q[1];
sx q[1];
rz(9.33334155975982) q[1];
cx q[1],q[0];
rz(-0.519732415676117) q[0];
sx q[0];
rz(3.87176528771455) q[0];
sx q[0];
rz(10.3768903970639) q[0];
rz(-2.04397630691528) q[2];
sx q[2];
rz(3.43103239138658) q[2];
sx q[2];
rz(12.0850665330808) q[2];
cx q[2],q[1];
rz(0.581272959709167) q[1];
sx q[1];
rz(2.88985073764855) q[1];
sx q[1];
rz(9.76141656040355) q[1];
rz(-0.594776928424835) q[3];
sx q[3];
rz(4.89434161980683) q[3];
sx q[3];
rz(9.86905491947337) q[3];
cx q[3],q[2];
rz(-0.664491951465607) q[2];
sx q[2];
rz(4.11509338219697) q[2];
sx q[2];
rz(10.5507803916852) q[2];
rz(-0.275157183408737) q[3];
sx q[3];
rz(3.75188437302644) q[3];
sx q[3];
rz(10.3278600931089) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.00984109565615654) q[0];
sx q[0];
rz(3.83509507973725) q[0];
sx q[0];
rz(10.118349826328) q[0];
rz(1.09610784053802) q[1];
sx q[1];
rz(4.1254417022043) q[1];
sx q[1];
rz(9.23446137308284) q[1];
cx q[1],q[0];
rz(0.411025583744049) q[0];
sx q[0];
rz(4.30982223351533) q[0];
sx q[0];
rz(10.7797540187757) q[0];
rz(0.960360407829285) q[2];
sx q[2];
rz(2.51548406680162) q[2];
sx q[2];
rz(9.92430630921527) q[2];
cx q[2],q[1];
rz(-0.438605010509491) q[1];
sx q[1];
rz(3.86400571663911) q[1];
sx q[1];
rz(9.45850678383514) q[1];
rz(1.05035972595215) q[3];
sx q[3];
rz(3.94030240376527) q[3];
sx q[3];
rz(9.91375899910136) q[3];
cx q[3],q[2];
rz(0.507431209087372) q[2];
sx q[2];
rz(3.85334417422349) q[2];
sx q[2];
rz(10.6445315837781) q[2];
rz(-0.142734199762344) q[3];
sx q[3];
rz(4.1536490042978) q[3];
sx q[3];
rz(10.3334096431653) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.742935657501221) q[0];
sx q[0];
rz(4.23149541218812) q[0];
sx q[0];
rz(10.2450364589612) q[0];
rz(-0.292076855897903) q[1];
sx q[1];
rz(4.21560362179811) q[1];
sx q[1];
rz(10.6727978944699) q[1];
cx q[1],q[0];
rz(0.672786355018616) q[0];
sx q[0];
rz(3.74816301663453) q[0];
sx q[0];
rz(10.8326921224515) q[0];
rz(1.36356472969055) q[2];
sx q[2];
rz(4.77522209485108) q[2];
sx q[2];
rz(9.83394763468906) q[2];
cx q[2],q[1];
rz(0.655908107757568) q[1];
sx q[1];
rz(3.81201330025727) q[1];
sx q[1];
rz(10.8877250909726) q[1];
rz(1.27078378200531) q[3];
sx q[3];
rz(4.08327016432817) q[3];
sx q[3];
rz(9.99216357468768) q[3];
cx q[3],q[2];
rz(-0.320392072200775) q[2];
sx q[2];
rz(4.19207874138887) q[2];
sx q[2];
rz(10.6382043123166) q[2];
rz(0.164725661277771) q[3];
sx q[3];
rz(4.03562596638734) q[3];
sx q[3];
rz(9.50638863294526) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.409590125083923) q[0];
sx q[0];
rz(4.42326834996278) q[0];
sx q[0];
rz(9.23871453701659) q[0];
rz(-0.204467251896858) q[1];
sx q[1];
rz(3.61354395945603) q[1];
sx q[1];
rz(11.2692334413449) q[1];
cx q[1],q[0];
rz(0.0649177953600883) q[0];
sx q[0];
rz(3.27900313039357) q[0];
sx q[0];
rz(9.74641836284801) q[0];
rz(0.517192780971527) q[2];
sx q[2];
rz(3.471374334889) q[2];
sx q[2];
rz(11.0512327909391) q[2];
cx q[2],q[1];
rz(-0.907210052013397) q[1];
sx q[1];
rz(2.13005903561647) q[1];
sx q[1];
rz(10.1480782389562) q[1];
rz(0.0760907977819443) q[3];
sx q[3];
rz(3.32167937059934) q[3];
sx q[3];
rz(10.5125125408094) q[3];
cx q[3],q[2];
rz(2.23564481735229) q[2];
sx q[2];
rz(3.93418953021104) q[2];
sx q[2];
rz(8.50558445452853) q[2];
rz(0.321331888437271) q[3];
sx q[3];
rz(4.20980957348878) q[3];
sx q[3];
rz(10.6725938081662) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.17386756837368) q[0];
sx q[0];
rz(3.64810213645036) q[0];
sx q[0];
rz(10.3029733061711) q[0];
rz(-1.3252660036087) q[1];
sx q[1];
rz(4.5538422187143) q[1];
sx q[1];
rz(10.8510659694593) q[1];
cx q[1],q[0];
rz(0.257098972797394) q[0];
sx q[0];
rz(4.63031938870484) q[0];
sx q[0];
rz(10.0819887876432) q[0];
rz(0.312517821788788) q[2];
sx q[2];
rz(4.48582819302613) q[2];
sx q[2];
rz(9.08337715863391) q[2];
cx q[2],q[1];
rz(0.519022762775421) q[1];
sx q[1];
rz(1.67913237412507) q[1];
sx q[1];
rz(10.1575570464055) q[1];
rz(1.24184668064117) q[3];
sx q[3];
rz(3.94021669228608) q[3];
sx q[3];
rz(9.19969520568057) q[3];
cx q[3],q[2];
rz(1.272012591362) q[2];
sx q[2];
rz(15*pi/13) q[2];
sx q[2];
rz(9.38298323600694) q[2];
rz(0.0614918656647205) q[3];
sx q[3];
rz(4.23065176804597) q[3];
sx q[3];
rz(9.86154705881282) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.35492658615112) q[0];
sx q[0];
rz(3.2275455092364) q[0];
sx q[0];
rz(10.45527110099) q[0];
rz(0.739739179611206) q[1];
sx q[1];
rz(4.67023530800874) q[1];
sx q[1];
rz(9.99634430407687) q[1];
cx q[1],q[0];
rz(0.334869027137756) q[0];
sx q[0];
rz(3.95831516583497) q[0];
sx q[0];
rz(10.263471043102) q[0];
rz(0.665847778320313) q[2];
sx q[2];
rz(4.81253496010835) q[2];
sx q[2];
rz(9.38728475048348) q[2];
cx q[2],q[1];
rz(0.953227698802948) q[1];
sx q[1];
rz(4.25848415692384) q[1];
sx q[1];
rz(10.7277165412824) q[1];
rz(0.86963301897049) q[3];
sx q[3];
rz(4.0881086905771) q[3];
sx q[3];
rz(9.55725505053207) q[3];
cx q[3],q[2];
rz(0.173438653349876) q[2];
sx q[2];
rz(3.88351586659486) q[2];
sx q[2];
rz(9.98897734879657) q[2];
rz(-0.126002222299576) q[3];
sx q[3];
rz(4.59994021256501) q[3];
sx q[3];
rz(9.71939394473239) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.09037041664124) q[0];
sx q[0];
rz(2.94169609447057) q[0];
sx q[0];
rz(10.1370545387189) q[0];
rz(-0.525871098041534) q[1];
sx q[1];
rz(3.55786782701547) q[1];
sx q[1];
rz(10.0902884364049) q[1];
cx q[1],q[0];
rz(0.804888010025024) q[0];
sx q[0];
rz(3.38382001419599) q[0];
sx q[0];
rz(9.42047591823294) q[0];
rz(0.697360694408417) q[2];
sx q[2];
rz(4.82657411892945) q[2];
sx q[2];
rz(9.6402028709571) q[2];
cx q[2],q[1];
rz(0.293840885162354) q[1];
sx q[1];
rz(4.49331513245637) q[1];
sx q[1];
rz(9.89311990737125) q[1];
rz(0.532879590988159) q[3];
sx q[3];
rz(2.86231923301751) q[3];
sx q[3];
rz(8.90002599953815) q[3];
cx q[3],q[2];
rz(0.157263815402985) q[2];
sx q[2];
rz(3.80165073473985) q[2];
sx q[2];
rz(10.8475793361585) q[2];
rz(-0.115193650126457) q[3];
sx q[3];
rz(3.65677496989305) q[3];
sx q[3];
rz(9.6173681974332) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.0499062351882458) q[0];
sx q[0];
rz(4.27727833588655) q[0];
sx q[0];
rz(9.23394749163791) q[0];
rz(0.626753687858582) q[1];
sx q[1];
rz(4.15767935116822) q[1];
sx q[1];
rz(9.7634886264722) q[1];
cx q[1],q[0];
rz(-0.209137976169586) q[0];
sx q[0];
rz(4.86281672318513) q[0];
sx q[0];
rz(9.77229533194705) q[0];
rz(2.023188829422) q[2];
sx q[2];
rz(3.63651108940179) q[2];
sx q[2];
rz(9.95646629332706) q[2];
cx q[2],q[1];
rz(-0.457044929265976) q[1];
sx q[1];
rz(4.23999741871888) q[1];
sx q[1];
rz(10.870949125282) q[1];
rz(2*pi/9) q[3];
sx q[3];
rz(3.55814576347406) q[3];
sx q[3];
rz(10.4884677886884) q[3];
cx q[3],q[2];
rz(0.646158576011658) q[2];
sx q[2];
rz(3.77604529460008) q[2];
sx q[2];
rz(8.72434190510913) q[2];
rz(0.897984802722931) q[3];
sx q[3];
rz(4.42823365529115) q[3];
sx q[3];
rz(9.59829100071594) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.531964659690857) q[0];
sx q[0];
rz(3.62319237192208) q[0];
sx q[0];
rz(10.0831260442655) q[0];
rz(0.610939621925354) q[1];
sx q[1];
rz(4.39716497262055) q[1];
sx q[1];
rz(9.56437408029243) q[1];
cx q[1],q[0];
rz(1.08482074737549) q[0];
sx q[0];
rz(3.07731008728082) q[0];
sx q[0];
rz(9.74150014518901) q[0];
rz(-1.62956500053406) q[2];
sx q[2];
rz(3.95707968075807) q[2];
sx q[2];
rz(11.3223222255628) q[2];
cx q[2],q[1];
rz(0.823972046375275) q[1];
sx q[1];
rz(3.66756996710832) q[1];
sx q[1];
rz(7.72061524390384) q[1];
rz(-0.0698088780045509) q[3];
sx q[3];
rz(3.78725615342195) q[3];
sx q[3];
rz(10.496746993057) q[3];
cx q[3],q[2];
rz(1.52479410171509) q[2];
sx q[2];
rz(4.15235880215699) q[2];
sx q[2];
rz(9.09366334079906) q[2];
rz(0.757742762565613) q[3];
sx q[3];
rz(3.53042200406129) q[3];
sx q[3];
rz(9.51265707462236) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.29632937908173) q[0];
sx q[0];
rz(4.15265754063661) q[0];
sx q[0];
rz(9.61035799085304) q[0];
rz(1.09620761871338) q[1];
sx q[1];
rz(3.35621598561341) q[1];
sx q[1];
rz(7.94013354777499) q[1];
cx q[1],q[0];
rz(2.05087041854858) q[0];
sx q[0];
rz(4.44286063511903) q[0];
sx q[0];
rz(8.96804497241184) q[0];
rz(0.276167780160904) q[2];
sx q[2];
rz(4.49509266217286) q[2];
sx q[2];
rz(10.9911918401639) q[2];
cx q[2],q[1];
rz(0.465860426425934) q[1];
sx q[1];
rz(4.61301127274568) q[1];
sx q[1];
rz(9.04870266317531) q[1];
rz(0.00580711988732219) q[3];
sx q[3];
rz(2.43666371901567) q[3];
sx q[3];
rz(10.1568069815557) q[3];
cx q[3],q[2];
rz(0.934027016162872) q[2];
sx q[2];
rz(4.03105387290055) q[2];
sx q[2];
rz(8.87256947755023) q[2];
rz(0.777837157249451) q[3];
sx q[3];
rz(3.99929556448991) q[3];
sx q[3];
rz(9.94267765282794) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.6388543844223) q[0];
sx q[0];
rz(3.32028447289998) q[0];
sx q[0];
rz(9.08799210786029) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(1.55953001976013) q[1];
sx q[1];
rz(3.96013704140718) q[1];
sx q[1];
rz(9.71163112520381) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(0.393331408500671) q[2];
sx q[2];
rz(4.03795102437074) q[2];
sx q[2];
rz(9.24436668156787) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(-0.460328578948975) q[3];
sx q[3];
rz(3.7478091438585) q[3];
sx q[3];
rz(10.6408948659818) q[3];
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
