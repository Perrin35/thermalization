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
rz(2.47406482696533) q[0];
sx q[0];
rz(2.3025371154123) q[0];
sx q[0];
rz(8.99230579137012) q[0];
rz(0.679838120937347) q[1];
sx q[1];
rz(4.82305839856202) q[1];
sx q[1];
rz(11.8001773118894) q[1];
cx q[1],q[0];
rz(1.95694088935852) q[0];
sx q[0];
rz(2.7650524695688) q[0];
sx q[0];
rz(9.17961054145499) q[0];
rz(-1.60176599025726) q[2];
sx q[2];
rz(8.03494325478608) q[2];
sx q[2];
rz(9.0860853254716) q[2];
cx q[2],q[1];
rz(2.74424648284912) q[1];
sx q[1];
rz(1.04342165787751) q[1];
sx q[1];
rz(11.1539476871411) q[1];
rz(7.00892639160156) q[3];
sx q[3];
rz(2.95765060384805) q[3];
sx q[3];
rz(8.20661268233463) q[3];
cx q[3],q[2];
rz(5.72797012329102) q[2];
sx q[2];
rz(3.74633553822572) q[2];
sx q[2];
rz(8.7629307269971) q[2];
rz(-1.23431789875031) q[3];
sx q[3];
rz(4.72287062008912) q[3];
sx q[3];
rz(9.57730566560432) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.90416514873505) q[0];
sx q[0];
rz(10.0769664367014) q[0];
sx q[0];
rz(12.1568393468778) q[0];
rz(2.88305902481079) q[1];
sx q[1];
rz(7.49001708825166) q[1];
sx q[1];
rz(9.97483394145175) q[1];
cx q[1],q[0];
rz(3.20018482208252) q[0];
sx q[0];
rz(0.482309015589305) q[0];
sx q[0];
rz(9.61260252296134) q[0];
rz(2.04437351226807) q[2];
sx q[2];
rz(5.22482219536836) q[2];
sx q[2];
rz(11.6771454572599) q[2];
cx q[2],q[1];
rz(0.0667958632111549) q[1];
sx q[1];
rz(2.65357119043405) q[1];
sx q[1];
rz(8.97843108176395) q[1];
rz(-0.370341777801514) q[3];
sx q[3];
rz(2.35622880061204) q[3];
sx q[3];
rz(11.5167431592862) q[3];
cx q[3],q[2];
rz(-0.516215801239014) q[2];
sx q[2];
rz(4.52793279488618) q[2];
sx q[2];
rz(11.974726176254) q[2];
rz(2.62434649467468) q[3];
sx q[3];
rz(5.17911687691743) q[3];
sx q[3];
rz(6.85784623622104) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(4.36106204986572) q[0];
sx q[0];
rz(4.99014201958711) q[0];
sx q[0];
rz(8.42380855082675) q[0];
rz(-1.71320223808289) q[1];
sx q[1];
rz(3.90053823788697) q[1];
sx q[1];
rz(15.7519745588224) q[1];
cx q[1],q[0];
rz(-1.50921702384949) q[0];
sx q[0];
rz(0.805018337565013) q[0];
sx q[0];
rz(10.4579143285672) q[0];
rz(0.736904561519623) q[2];
sx q[2];
rz(4.67469635804231) q[2];
sx q[2];
rz(8.4287615775983) q[2];
cx q[2],q[1];
rz(-0.219798669219017) q[1];
sx q[1];
rz(4.21854070027406) q[1];
sx q[1];
rz(7.06755683421298) q[1];
rz(-0.229960143566132) q[3];
sx q[3];
rz(0.176475675898143) q[3];
sx q[3];
rz(6.26173660754367) q[3];
cx q[3],q[2];
rz(-0.222378045320511) q[2];
sx q[2];
rz(4.69043758709962) q[2];
sx q[2];
rz(12.456557250015) q[2];
rz(-3.77171969413757) q[3];
sx q[3];
rz(2.58702418406541) q[3];
sx q[3];
rz(9.71234587430164) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.5146062374115) q[0];
sx q[0];
rz(2.08727482159669) q[0];
sx q[0];
rz(10.5188033342282) q[0];
rz(-0.473597228527069) q[1];
sx q[1];
rz(2.43515989382798) q[1];
sx q[1];
rz(11.8461341619413) q[1];
cx q[1],q[0];
rz(2.73137736320496) q[0];
sx q[0];
rz(3.71772989829118) q[0];
sx q[0];
rz(10.4124656677167) q[0];
rz(-4.85910654067993) q[2];
sx q[2];
rz(-0.311390487355641) q[2];
sx q[2];
rz(13.9420799970548) q[2];
cx q[2],q[1];
rz(-4.02948808670044) q[1];
sx q[1];
rz(1.27165368397767) q[1];
sx q[1];
rz(10.2350601315419) q[1];
rz(1.88915634155273) q[3];
sx q[3];
rz(5.44252458413179) q[3];
sx q[3];
rz(11.5454032182614) q[3];
cx q[3],q[2];
rz(-2.58064103126526) q[2];
sx q[2];
rz(1.66839423974092) q[2];
sx q[2];
rz(9.57693835198089) q[2];
rz(1.84693193435669) q[3];
sx q[3];
rz(4.84681645234162) q[3];
sx q[3];
rz(10.7919599771421) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(3.21402859687805) q[0];
sx q[0];
rz(2.57885882456834) q[0];
sx q[0];
rz(9.80612341164752) q[0];
rz(-0.588224768638611) q[1];
sx q[1];
rz(4.87600055535371) q[1];
sx q[1];
rz(9.3605819851081) q[1];
cx q[1],q[0];
rz(1.67468154430389) q[0];
sx q[0];
rz(4.88962522347505) q[0];
sx q[0];
rz(10.6195090770642) q[0];
rz(1.02673304080963) q[2];
sx q[2];
rz(3.83260366519029) q[2];
sx q[2];
rz(12.9071049451749) q[2];
cx q[2],q[1];
rz(-1.20926797389984) q[1];
sx q[1];
rz(1.28636828263337) q[1];
sx q[1];
rz(11.8032164335172) q[1];
rz(3.34745335578918) q[3];
sx q[3];
rz(1.69042173226411) q[3];
sx q[3];
rz(8.69900194405719) q[3];
cx q[3],q[2];
rz(-0.753078162670135) q[2];
sx q[2];
rz(4.45191028912599) q[2];
sx q[2];
rz(13.3591673135678) q[2];
rz(-0.134752839803696) q[3];
sx q[3];
rz(5.71767178376252) q[3];
sx q[3];
rz(7.77367148398563) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.438692688941956) q[0];
sx q[0];
rz(4.7287861426645) q[0];
sx q[0];
rz(8.47184274195834) q[0];
rz(1.89035248756409) q[1];
sx q[1];
rz(4.82149175007875) q[1];
sx q[1];
rz(6.04891893862888) q[1];
cx q[1],q[0];
rz(-1.33279848098755) q[0];
sx q[0];
rz(7.9991981108957) q[0];
sx q[0];
rz(10.4893735408704) q[0];
rz(-2.29069828987122) q[2];
sx q[2];
rz(5.70514455636079) q[2];
sx q[2];
rz(11.243779039375) q[2];
cx q[2],q[1];
rz(-1.79978787899017) q[1];
sx q[1];
rz(4.82625666459138) q[1];
sx q[1];
rz(13.1525428056638) q[1];
rz(-3.02604699134827) q[3];
sx q[3];
rz(6.92096391518647) q[3];
sx q[3];
rz(7.7496136188428) q[3];
cx q[3],q[2];
rz(0.171440631151199) q[2];
sx q[2];
rz(2.54153874714906) q[2];
sx q[2];
rz(8.47792652844592) q[2];
rz(-1.44174492359161) q[3];
sx q[3];
rz(4.46906057198579) q[3];
sx q[3];
rz(9.54615548103257) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.22987188398838) q[0];
sx q[0];
rz(2.13891020615632) q[0];
sx q[0];
rz(6.65807840823337) q[0];
rz(-3.88453841209412) q[1];
sx q[1];
rz(3.96447059710557) q[1];
sx q[1];
rz(8.76854149102374) q[1];
cx q[1],q[0];
rz(0.292268812656403) q[0];
sx q[0];
rz(4.97172299225862) q[0];
sx q[0];
rz(10.3056429982106) q[0];
rz(-4.07915449142456) q[2];
sx q[2];
rz(2.04549625714357) q[2];
sx q[2];
rz(8.43881717919513) q[2];
cx q[2],q[1];
rz(-3.37486839294434) q[1];
sx q[1];
rz(5.08298757870729) q[1];
sx q[1];
rz(7.40335772036716) q[1];
rz(-0.979009568691254) q[3];
sx q[3];
rz(4.87866345246369) q[3];
sx q[3];
rz(8.72445068358585) q[3];
cx q[3],q[2];
rz(-0.528342127799988) q[2];
sx q[2];
rz(0.873177917795726) q[2];
sx q[2];
rz(13.2163216829221) q[2];
rz(0.849082708358765) q[3];
sx q[3];
rz(5.64571061928804) q[3];
sx q[3];
rz(8.23684833048984) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.77173018455505) q[0];
sx q[0];
rz(4.71314612229402) q[0];
sx q[0];
rz(9.74298096298381) q[0];
rz(0.33266207575798) q[1];
sx q[1];
rz(3.6948488672548) q[1];
sx q[1];
rz(6.66019318103) q[1];
cx q[1],q[0];
rz(0.724125981330872) q[0];
sx q[0];
rz(1.23403492768342) q[0];
sx q[0];
rz(11.0444321393888) q[0];
rz(2.67443561553955) q[2];
sx q[2];
rz(4.74434390862519) q[2];
sx q[2];
rz(11.6866452455442) q[2];
cx q[2],q[1];
rz(0.748063802719116) q[1];
sx q[1];
rz(8.22581973870332) q[1];
sx q[1];
rz(10.2753690838735) q[1];
rz(-2.92253232002258) q[3];
sx q[3];
rz(5.78997603257234) q[3];
sx q[3];
rz(8.40531704425021) q[3];
cx q[3],q[2];
rz(0.134100213646889) q[2];
sx q[2];
rz(5.00517311890657) q[2];
sx q[2];
rz(7.21416399478122) q[2];
rz(-1.53829860687256) q[3];
sx q[3];
rz(1.80360165436799) q[3];
sx q[3];
rz(4.45796821116611) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.613049209117889) q[0];
sx q[0];
rz(0.407540233927318) q[0];
sx q[0];
rz(13.0079324006955) q[0];
rz(3.40670561790466) q[1];
sx q[1];
rz(4.84129896958406) q[1];
sx q[1];
rz(10.9354839086454) q[1];
cx q[1],q[0];
rz(-1.85454523563385) q[0];
sx q[0];
rz(2.00187388260896) q[0];
sx q[0];
rz(12.8104684114377) q[0];
rz(3.73305082321167) q[2];
sx q[2];
rz(4.42902127106721) q[2];
sx q[2];
rz(10.7705332994382) q[2];
cx q[2],q[1];
rz(5.93055486679077) q[1];
sx q[1];
rz(4.81204679806764) q[1];
sx q[1];
rz(9.65865512787505) q[1];
rz(-3.61169958114624) q[3];
sx q[3];
rz(10.130887182551) q[3];
sx q[3];
rz(11.6942672490995) q[3];
cx q[3],q[2];
rz(-0.604921698570251) q[2];
sx q[2];
rz(4.41854718525941) q[2];
sx q[2];
rz(9.89337975382015) q[2];
rz(2.12971496582031) q[3];
sx q[3];
rz(4.56465819676454) q[3];
sx q[3];
rz(8.34356305598422) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.68779468536377) q[0];
sx q[0];
rz(3.98465487559373) q[0];
sx q[0];
rz(11.1296144485395) q[0];
rz(-0.948024809360504) q[1];
sx q[1];
rz(4.94322315056855) q[1];
sx q[1];
rz(12.2925913095395) q[1];
cx q[1],q[0];
rz(-1.70234668254852) q[0];
sx q[0];
rz(2.99860234757001) q[0];
sx q[0];
rz(8.37583491801425) q[0];
rz(4.06616973876953) q[2];
sx q[2];
rz(4.71018937428529) q[2];
sx q[2];
rz(8.45797214507266) q[2];
cx q[2],q[1];
rz(5.93310165405273) q[1];
sx q[1];
rz(5.5701688845926) q[1];
sx q[1];
rz(16.5258364438932) q[1];
rz(0.540739774703979) q[3];
sx q[3];
rz(2.39168366988237) q[3];
sx q[3];
rz(12.03394982814) q[3];
cx q[3],q[2];
rz(0.561293661594391) q[2];
sx q[2];
rz(4.21458474000032) q[2];
sx q[2];
rz(13.0143558740537) q[2];
rz(2.50428462028503) q[3];
sx q[3];
rz(4.16770485241944) q[3];
sx q[3];
rz(8.52972856759235) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.238272130489349) q[0];
sx q[0];
rz(0.11549177964265) q[0];
sx q[0];
rz(9.83066312073871) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(4.38971662521362) q[1];
sx q[1];
rz(4.53563955624635) q[1];
sx q[1];
rz(5.44967005252048) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(-5.85205221176147) q[2];
sx q[2];
rz(1.7196503003412) q[2];
sx q[2];
rz(15.7859973668973) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(4.29893827438354) q[3];
sx q[3];
rz(5.754383238154) q[3];
sx q[3];
rz(12.567440724365) q[3];
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
