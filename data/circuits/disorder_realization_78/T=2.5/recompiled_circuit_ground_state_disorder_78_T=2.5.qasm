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
rz(2.65277695655823) q[0];
sx q[0];
rz(4.10696944792802) q[0];
sx q[0];
rz(10.3416359186093) q[0];
rz(1.44044530391693) q[1];
sx q[1];
rz(2.11357656319673) q[1];
sx q[1];
rz(10.1205747485082) q[1];
cx q[1],q[0];
rz(-1.35499334335327) q[0];
sx q[0];
rz(3.57172122796113) q[0];
sx q[0];
rz(9.29456651806041) q[0];
rz(-1.53343391418457) q[2];
sx q[2];
rz(4.76000216801698) q[2];
sx q[2];
rz(11.1936932563703) q[2];
cx q[2],q[1];
rz(4.02278089523315) q[1];
sx q[1];
rz(5.08663740952546) q[1];
sx q[1];
rz(9.24140518008872) q[1];
rz(2.51812434196472) q[3];
sx q[3];
rz(0.586373241739818) q[3];
sx q[3];
rz(9.24345195888683) q[3];
cx q[3],q[2];
rz(0.157713890075684) q[2];
sx q[2];
rz(2.32150510151918) q[2];
sx q[2];
rz(8.52436707019016) q[2];
rz(-0.793120980262756) q[3];
sx q[3];
rz(4.79149881203706) q[3];
sx q[3];
rz(10.1233842730443) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.52914905548096) q[0];
sx q[0];
rz(1.97310331662232) q[0];
sx q[0];
rz(10.1103041529576) q[0];
rz(2.70870518684387) q[1];
sx q[1];
rz(1.87501779397065) q[1];
sx q[1];
rz(8.11087045668765) q[1];
cx q[1],q[0];
rz(3.20852422714233) q[0];
sx q[0];
rz(4.21837225754792) q[0];
sx q[0];
rz(10.0426725506703) q[0];
rz(-0.493538022041321) q[2];
sx q[2];
rz(3.94935777981813) q[2];
sx q[2];
rz(8.27880630492374) q[2];
cx q[2],q[1];
rz(-0.594562768936157) q[1];
sx q[1];
rz(0.702077539759227) q[1];
sx q[1];
rz(11.8005676031034) q[1];
rz(-1.01382970809937) q[3];
sx q[3];
rz(1.7855705340677) q[3];
sx q[3];
rz(11.545629477493) q[3];
cx q[3],q[2];
rz(0.965852975845337) q[2];
sx q[2];
rz(2.73128351767594) q[2];
sx q[2];
rz(9.1773447304885) q[2];
rz(-0.115821726620197) q[3];
sx q[3];
rz(4.85627928574617) q[3];
sx q[3];
rz(10.003487443916) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.98473572731018) q[0];
sx q[0];
rz(5.90710154374177) q[0];
sx q[0];
rz(10.4896520137708) q[0];
rz(-0.633060097694397) q[1];
sx q[1];
rz(5.12448373635346) q[1];
sx q[1];
rz(12.4298248052518) q[1];
cx q[1],q[0];
rz(-0.318131923675537) q[0];
sx q[0];
rz(5.54654351075227) q[0];
sx q[0];
rz(7.89026854037448) q[0];
rz(0.630332410335541) q[2];
sx q[2];
rz(3.69519153435762) q[2];
sx q[2];
rz(11.8198616266172) q[2];
cx q[2],q[1];
rz(-0.284408390522003) q[1];
sx q[1];
rz(2.29400024016435) q[1];
sx q[1];
rz(13.4884233236234) q[1];
rz(-0.430604159832001) q[3];
sx q[3];
rz(4.33093682129914) q[3];
sx q[3];
rz(9.8051999270837) q[3];
cx q[3],q[2];
rz(2.44240498542786) q[2];
sx q[2];
rz(5.03085378010804) q[2];
sx q[2];
rz(7.88318071364566) q[2];
rz(-1.09347462654114) q[3];
sx q[3];
rz(1.63960614998872) q[3];
sx q[3];
rz(9.77017412184879) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.714447557926178) q[0];
sx q[0];
rz(3.40304279525811) q[0];
sx q[0];
rz(10.0015883803289) q[0];
rz(3.86911129951477) q[1];
sx q[1];
rz(5.19176879723603) q[1];
sx q[1];
rz(11.6980325937192) q[1];
cx q[1],q[0];
rz(1.51423466205597) q[0];
sx q[0];
rz(4.31656506856019) q[0];
sx q[0];
rz(8.53781441449329) q[0];
rz(2.5593535900116) q[2];
sx q[2];
rz(4.2622196992212) q[2];
sx q[2];
rz(5.76046392916843) q[2];
cx q[2],q[1];
rz(-2.78177452087402) q[1];
sx q[1];
rz(5.0553698857599) q[1];
sx q[1];
rz(10.305007493488) q[1];
rz(1.06290137767792) q[3];
sx q[3];
rz(1.00818553765351) q[3];
sx q[3];
rz(12.3131761312406) q[3];
cx q[3],q[2];
rz(3.33680248260498) q[2];
sx q[2];
rz(1.75274399121339) q[2];
sx q[2];
rz(8.72926572560474) q[2];
rz(2.92693114280701) q[3];
sx q[3];
rz(4.64273944695527) q[3];
sx q[3];
rz(13.3796794176023) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.60140860080719) q[0];
sx q[0];
rz(3.75647255976731) q[0];
sx q[0];
rz(11.4662577867429) q[0];
rz(2.93101382255554) q[1];
sx q[1];
rz(5.6497804244333) q[1];
sx q[1];
rz(11.0731552600782) q[1];
cx q[1],q[0];
rz(-0.570151448249817) q[0];
sx q[0];
rz(4.23415830929811) q[0];
sx q[0];
rz(11.5290407895963) q[0];
rz(2.9635648727417) q[2];
sx q[2];
rz(4.79498663743074) q[2];
sx q[2];
rz(11.6043932199399) q[2];
cx q[2],q[1];
rz(2.32012176513672) q[1];
sx q[1];
rz(4.55685249169404) q[1];
sx q[1];
rz(8.93225896953746) q[1];
rz(4.04519653320313) q[3];
sx q[3];
rz(1.76931873162324) q[3];
sx q[3];
rz(9.0302607178609) q[3];
cx q[3],q[2];
rz(1.41399991512299) q[2];
sx q[2];
rz(1.71213451226289) q[2];
sx q[2];
rz(8.52777907847568) q[2];
rz(-1.2796847820282) q[3];
sx q[3];
rz(1.01062670548493) q[3];
sx q[3];
rz(9.47198734655186) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.14192640781403) q[0];
sx q[0];
rz(3.76745578845079) q[0];
sx q[0];
rz(9.75489193796321) q[0];
rz(5.64239549636841) q[1];
sx q[1];
rz(4.90744129021699) q[1];
sx q[1];
rz(11.5908367395322) q[1];
cx q[1],q[0];
rz(0.278110891580582) q[0];
sx q[0];
rz(4.76910284359986) q[0];
sx q[0];
rz(11.5866505861203) q[0];
rz(2.44637441635132) q[2];
sx q[2];
rz(0.406031282740184) q[2];
sx q[2];
rz(6.75008127688571) q[2];
cx q[2],q[1];
rz(4.35738039016724) q[1];
sx q[1];
rz(6.81095710595185) q[1];
sx q[1];
rz(13.5548262357633) q[1];
rz(0.401310741901398) q[3];
sx q[3];
rz(5.82914033730561) q[3];
sx q[3];
rz(11.2140257120053) q[3];
cx q[3],q[2];
rz(-1.28341245651245) q[2];
sx q[2];
rz(3.39571839769418) q[2];
sx q[2];
rz(8.53509429692432) q[2];
rz(2.54559183120728) q[3];
sx q[3];
rz(4.21079769928987) q[3];
sx q[3];
rz(9.53143552540942) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.135705307126045) q[0];
sx q[0];
rz(4.27771797974641) q[0];
sx q[0];
rz(8.29179046153232) q[0];
rz(0.639444470405579) q[1];
sx q[1];
rz(2.06381276448304) q[1];
sx q[1];
rz(11.5989451169889) q[1];
cx q[1],q[0];
rz(3.29776906967163) q[0];
sx q[0];
rz(1.42089095910127) q[0];
sx q[0];
rz(8.90971533059284) q[0];
rz(0.329154700040817) q[2];
sx q[2];
rz(5.07369306881959) q[2];
sx q[2];
rz(9.87701723574802) q[2];
cx q[2],q[1];
rz(2.1167368888855) q[1];
sx q[1];
rz(4.9369764645868) q[1];
sx q[1];
rz(9.54292407482072) q[1];
rz(-0.308150917291641) q[3];
sx q[3];
rz(2.2053092439943) q[3];
sx q[3];
rz(9.24361678063079) q[3];
cx q[3],q[2];
rz(1.39666676521301) q[2];
sx q[2];
rz(5.27417651017243) q[2];
sx q[2];
rz(10.3691150307576) q[2];
rz(0.622632920742035) q[3];
sx q[3];
rz(4.12992009718949) q[3];
sx q[3];
rz(8.6377641916196) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(3.0866162776947) q[0];
sx q[0];
rz(3.94785461028154) q[0];
sx q[0];
rz(9.42900931266277) q[0];
rz(1.46842336654663) q[1];
sx q[1];
rz(5.51989689667756) q[1];
sx q[1];
rz(9.25031459926769) q[1];
cx q[1],q[0];
rz(-2.3136351108551) q[0];
sx q[0];
rz(3.35010977287824) q[0];
sx q[0];
rz(9.88693342208072) q[0];
rz(-1.78082990646362) q[2];
sx q[2];
rz(5.82593432267243) q[2];
sx q[2];
rz(8.78733251094028) q[2];
cx q[2],q[1];
rz(1.46943879127502) q[1];
sx q[1];
rz(4.64635196526582) q[1];
sx q[1];
rz(10.1237453579824) q[1];
rz(0.264411926269531) q[3];
sx q[3];
rz(3.38591915567453) q[3];
sx q[3];
rz(8.64230409859821) q[3];
cx q[3],q[2];
rz(-1.89094948768616) q[2];
sx q[2];
rz(6.18929425080354) q[2];
sx q[2];
rz(8.88442144393131) q[2];
rz(-0.0171026084572077) q[3];
sx q[3];
rz(0.981868656473704) q[3];
sx q[3];
rz(8.76300022601291) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.09062016010284) q[0];
sx q[0];
rz(7.08565917809541) q[0];
sx q[0];
rz(10.2113766431729) q[0];
rz(-2.03589177131653) q[1];
sx q[1];
rz(1.66037550766999) q[1];
sx q[1];
rz(12.3451442480008) q[1];
cx q[1],q[0];
rz(1.82411372661591) q[0];
sx q[0];
rz(5.41659441788728) q[0];
sx q[0];
rz(11.3636092901151) q[0];
rz(-3.31516289710999) q[2];
sx q[2];
rz(4.61887958844239) q[2];
sx q[2];
rz(11.0509183168332) q[2];
cx q[2],q[1];
rz(2.55334997177124) q[1];
sx q[1];
rz(2.60946312745149) q[1];
sx q[1];
rz(8.27178070544406) q[1];
rz(0.880582988262177) q[3];
sx q[3];
rz(4.89403787453706) q[3];
sx q[3];
rz(12.0491373300473) q[3];
cx q[3],q[2];
rz(-1.14563870429993) q[2];
sx q[2];
rz(6.52572408516938) q[2];
sx q[2];
rz(7.46875927447482) q[2];
rz(2.03541851043701) q[3];
sx q[3];
rz(5.26148024399812) q[3];
sx q[3];
rz(8.4405881524007) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.361586540937424) q[0];
sx q[0];
rz(4.92585697968537) q[0];
sx q[0];
rz(10.1345661640088) q[0];
rz(-2.24787259101868) q[1];
sx q[1];
rz(3.76944670279557) q[1];
sx q[1];
rz(9.89178017377063) q[1];
cx q[1],q[0];
rz(-0.929995357990265) q[0];
sx q[0];
rz(3.44307160575921) q[0];
sx q[0];
rz(8.94433904289409) q[0];
rz(-1.30206656455994) q[2];
sx q[2];
rz(4.10435697634751) q[2];
sx q[2];
rz(15.7096085309903) q[2];
cx q[2],q[1];
rz(-2.26450061798096) q[1];
sx q[1];
rz(2.90007259150083) q[1];
sx q[1];
rz(15.5515222310941) q[1];
rz(0.696749627590179) q[3];
sx q[3];
rz(1.02068868477876) q[3];
sx q[3];
rz(10.1969800949018) q[3];
cx q[3],q[2];
rz(-1.34221565723419) q[2];
sx q[2];
rz(4.04664692481095) q[2];
sx q[2];
rz(9.56077710389301) q[2];
rz(0.329867660999298) q[3];
sx q[3];
rz(2.07432750065858) q[3];
sx q[3];
rz(9.72190386652156) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-2.57294774055481) q[0];
sx q[0];
rz(5.26684180100495) q[0];
sx q[0];
rz(9.96799216269656) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(-2.41653227806091) q[1];
sx q[1];
rz(1.44581702549989) q[1];
sx q[1];
rz(12.8952305078427) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(-1.8419908285141) q[2];
sx q[2];
rz(3.35077236791188) q[2];
sx q[2];
rz(10.8727918624799) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(1.60652339458466) q[3];
sx q[3];
rz(2.01051452954347) q[3];
sx q[3];
rz(10.1955541729848) q[3];
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
