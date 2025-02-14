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
rz(-0.152205169200897) q[0];
sx q[0];
rz(2.92899911304051) q[0];
sx q[0];
rz(8.69060353039905) q[0];
rz(-1.92501592636108) q[1];
sx q[1];
rz(2.79563784797723) q[1];
sx q[1];
rz(12.5911271333615) q[1];
cx q[1],q[0];
rz(0.615413248538971) q[0];
sx q[0];
rz(3.21747412730987) q[0];
sx q[0];
rz(11.6021468400876) q[0];
rz(3.58777093887329) q[2];
sx q[2];
rz(6.02503743966157) q[2];
sx q[2];
rz(6.4050810098569) q[2];
cx q[2],q[1];
rz(-0.369444459676743) q[1];
sx q[1];
rz(1.27190628846223) q[1];
sx q[1];
rz(12.6211838483731) q[1];
rz(0.253938943147659) q[3];
sx q[3];
rz(3.36081664462621) q[3];
sx q[3];
rz(8.49646470545932) q[3];
cx q[3],q[2];
rz(-0.884322464466095) q[2];
sx q[2];
rz(1.49296441872651) q[2];
sx q[2];
rz(12.2575447320859) q[2];
rz(-2.31579995155334) q[3];
sx q[3];
rz(4.14958575566346) q[3];
sx q[3];
rz(8.66081658600971) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.170771330595016) q[0];
sx q[0];
rz(5.05404833157594) q[0];
sx q[0];
rz(7.31087038516208) q[0];
rz(0.816142559051514) q[1];
sx q[1];
rz(7.40152421792085) q[1];
sx q[1];
rz(8.61229387520953) q[1];
cx q[1],q[0];
rz(3.66175627708435) q[0];
sx q[0];
rz(4.87126377423341) q[0];
sx q[0];
rz(12.2126416921537) q[0];
rz(1.19742894172668) q[2];
sx q[2];
rz(4.72561422188813) q[2];
sx q[2];
rz(11.6475078821103) q[2];
cx q[2],q[1];
rz(-2.83904385566711) q[1];
sx q[1];
rz(4.85022071202333) q[1];
sx q[1];
rz(12.0518100023191) q[1];
rz(-1.22187542915344) q[3];
sx q[3];
rz(8.18980613549287) q[3];
sx q[3];
rz(9.27415098845168) q[3];
cx q[3],q[2];
rz(0.574277877807617) q[2];
sx q[2];
rz(1.48941663106019) q[2];
sx q[2];
rz(8.50585023163959) q[2];
rz(0.531416654586792) q[3];
sx q[3];
rz(1.04549494584138) q[3];
sx q[3];
rz(9.46596379055783) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.703192174434662) q[0];
sx q[0];
rz(4.21125355561311) q[0];
sx q[0];
rz(10.5620919227521) q[0];
rz(4.89262676239014) q[1];
sx q[1];
rz(0.556823404627391) q[1];
sx q[1];
rz(10.1577387213628) q[1];
cx q[1],q[0];
rz(-1.83702862262726) q[0];
sx q[0];
rz(3.65849980910356) q[0];
sx q[0];
rz(7.23148629664584) q[0];
rz(-0.60422271490097) q[2];
sx q[2];
rz(5.19622031052644) q[2];
sx q[2];
rz(8.41755816935703) q[2];
cx q[2],q[1];
rz(0.939483523368835) q[1];
sx q[1];
rz(2.04014781315858) q[1];
sx q[1];
rz(8.34581158160373) q[1];
rz(2.69286608695984) q[3];
sx q[3];
rz(4.70646193821961) q[3];
sx q[3];
rz(11.4420044183652) q[3];
cx q[3],q[2];
rz(-0.00400323374196887) q[2];
sx q[2];
rz(4.87103775342042) q[2];
sx q[2];
rz(8.3859046459119) q[2];
rz(2.13933563232422) q[3];
sx q[3];
rz(4.44332364399964) q[3];
sx q[3];
rz(5.99302241801425) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.206594735383987) q[0];
sx q[0];
rz(4.22797122796113) q[0];
sx q[0];
rz(8.51759061812564) q[0];
rz(-0.157792374491692) q[1];
sx q[1];
rz(2.12674990494783) q[1];
sx q[1];
rz(8.14354572295352) q[1];
cx q[1],q[0];
rz(0.322037547826767) q[0];
sx q[0];
rz(1.05196181138093) q[0];
sx q[0];
rz(10.2793209314267) q[0];
rz(-1.96958136558533) q[2];
sx q[2];
rz(5.38790169556672) q[2];
sx q[2];
rz(9.14384732245609) q[2];
cx q[2],q[1];
rz(-2.70420265197754) q[1];
sx q[1];
rz(10.5713974555307) q[1];
sx q[1];
rz(8.54027042388126) q[1];
rz(-5.53604936599731) q[3];
sx q[3];
rz(7.1077562888437) q[3];
sx q[3];
rz(7.41546747683688) q[3];
cx q[3],q[2];
rz(-0.247187569737434) q[2];
sx q[2];
rz(2.20535680850083) q[2];
sx q[2];
rz(10.7144942045133) q[2];
rz(-1.79734086990356) q[3];
sx q[3];
rz(4.59857633908326) q[3];
sx q[3];
rz(10.124919986717) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.460713893175125) q[0];
sx q[0];
rz(2.64321556885774) q[0];
sx q[0];
rz(7.30145308970615) q[0];
rz(1.10300385951996) q[1];
sx q[1];
rz(3.85351255734498) q[1];
sx q[1];
rz(8.08432254790469) q[1];
cx q[1],q[0];
rz(-1.06062746047974) q[0];
sx q[0];
rz(0.617653520899363) q[0];
sx q[0];
rz(7.20213315486118) q[0];
rz(1.81183397769928) q[2];
sx q[2];
rz(2.75169810851152) q[2];
sx q[2];
rz(9.76788536309406) q[2];
cx q[2],q[1];
rz(-0.835070431232452) q[1];
sx q[1];
rz(4.84185782273347) q[1];
sx q[1];
rz(7.7259217262189) q[1];
rz(-0.540328979492188) q[3];
sx q[3];
rz(4.76183882554109) q[3];
sx q[3];
rz(10.845239853851) q[3];
cx q[3],q[2];
rz(1.54201424121857) q[2];
sx q[2];
rz(5.3946489413553) q[2];
sx q[2];
rz(7.28518769740268) q[2];
rz(0.691397845745087) q[3];
sx q[3];
rz(4.32203558285768) q[3];
sx q[3];
rz(10.9455538749616) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.602867603302) q[0];
sx q[0];
rz(2.0724916775995) q[0];
sx q[0];
rz(10.050050353996) q[0];
rz(1.94847762584686) q[1];
sx q[1];
rz(3.83139756520326) q[1];
sx q[1];
rz(8.85745737551852) q[1];
cx q[1],q[0];
rz(1.75174868106842) q[0];
sx q[0];
rz(4.53671828110749) q[0];
sx q[0];
rz(11.782920575134) q[0];
rz(3.88598966598511) q[2];
sx q[2];
rz(4.70687654812867) q[2];
sx q[2];
rz(6.65276930331394) q[2];
cx q[2],q[1];
rz(-0.159298270940781) q[1];
sx q[1];
rz(4.68641534646089) q[1];
sx q[1];
rz(10.0433071017186) q[1];
rz(-0.179320514202118) q[3];
sx q[3];
rz(8.02743545373017) q[3];
sx q[3];
rz(11.1389306545179) q[3];
cx q[3],q[2];
rz(0.83711963891983) q[2];
sx q[2];
rz(4.42507985432679) q[2];
sx q[2];
rz(4.72338960169956) q[2];
rz(2.52901530265808) q[3];
sx q[3];
rz(4.30227568944032) q[3];
sx q[3];
rz(9.26900560259029) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.09538578987122) q[0];
sx q[0];
rz(5.09236565430696) q[0];
sx q[0];
rz(10.5516315460126) q[0];
rz(4.41121053695679) q[1];
sx q[1];
rz(4.32952252228791) q[1];
sx q[1];
rz(9.94584009646579) q[1];
cx q[1],q[0];
rz(-1.23494863510132) q[0];
sx q[0];
rz(2.36813304026658) q[0];
sx q[0];
rz(9.93585882186099) q[0];
rz(-0.0959224551916122) q[2];
sx q[2];
rz(6.61357894738252) q[2];
sx q[2];
rz(10.4411972522657) q[2];
cx q[2],q[1];
rz(-1.48496294021606) q[1];
sx q[1];
rz(5.28826013405854) q[1];
sx q[1];
rz(13.5055718183438) q[1];
rz(-1.02977430820465) q[3];
sx q[3];
rz(4.47458270390565) q[3];
sx q[3];
rz(10.4600898981015) q[3];
cx q[3],q[2];
rz(2.11789393424988) q[2];
sx q[2];
rz(4.50569537480409) q[2];
sx q[2];
rz(10.6517395734708) q[2];
rz(1.46206939220428) q[3];
sx q[3];
rz(4.65897289116914) q[3];
sx q[3];
rz(9.86633447407886) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.113144971430302) q[0];
sx q[0];
rz(4.17130175431306) q[0];
sx q[0];
rz(8.87766215800449) q[0];
rz(0.0881349146366119) q[1];
sx q[1];
rz(4.93556764920289) q[1];
sx q[1];
rz(9.84731242655917) q[1];
cx q[1],q[0];
rz(-0.201326578855515) q[0];
sx q[0];
rz(2.8597587962919) q[0];
sx q[0];
rz(12.4875960111539) q[0];
rz(-0.369341522455215) q[2];
sx q[2];
rz(5.20526591141755) q[2];
sx q[2];
rz(11.7950238943021) q[2];
cx q[2],q[1];
rz(-0.402289807796478) q[1];
sx q[1];
rz(4.76843586762483) q[1];
sx q[1];
rz(9.35561898945972) q[1];
rz(0.651106774806976) q[3];
sx q[3];
rz(2.91251787741716) q[3];
sx q[3];
rz(9.56292911469146) q[3];
cx q[3],q[2];
rz(-1.22534584999084) q[2];
sx q[2];
rz(4.64606717427308) q[2];
sx q[2];
rz(9.71205151676341) q[2];
rz(0.542879939079285) q[3];
sx q[3];
rz(3.91176137526567) q[3];
sx q[3];
rz(9.36667536794349) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.74202370643616) q[0];
sx q[0];
rz(5.5990733226114) q[0];
sx q[0];
rz(10.2020887494008) q[0];
rz(2.21513915061951) q[1];
sx q[1];
rz(1.99942055543) q[1];
sx q[1];
rz(11.171411728851) q[1];
cx q[1],q[0];
rz(2.1374716758728) q[0];
sx q[0];
rz(3.65411940415437) q[0];
sx q[0];
rz(8.73449698685809) q[0];
rz(2.16996169090271) q[2];
sx q[2];
rz(4.92044320900971) q[2];
sx q[2];
rz(11.9171981573026) q[2];
cx q[2],q[1];
rz(1.39231264591217) q[1];
sx q[1];
rz(3.71807775099809) q[1];
sx q[1];
rz(6.91872904299899) q[1];
rz(2.18782210350037) q[3];
sx q[3];
rz(4.16954735119874) q[3];
sx q[3];
rz(11.4092464208524) q[3];
cx q[3],q[2];
rz(-0.425614386796951) q[2];
sx q[2];
rz(4.51407423813874) q[2];
sx q[2];
rz(7.37374899386569) q[2];
rz(0.9645094871521) q[3];
sx q[3];
rz(4.15300002892549) q[3];
sx q[3];
rz(8.20961210726901) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.575080931186676) q[0];
sx q[0];
rz(3.48001122673089) q[0];
sx q[0];
rz(10.9348498344342) q[0];
rz(0.0343200229108334) q[1];
sx q[1];
rz(1.80205539067323) q[1];
sx q[1];
rz(6.71520302294894) q[1];
cx q[1],q[0];
rz(0.851381897926331) q[0];
sx q[0];
rz(3.68358931143815) q[0];
sx q[0];
rz(9.53882286547824) q[0];
rz(2.92145419120789) q[2];
sx q[2];
rz(5.92414489586885) q[2];
sx q[2];
rz(9.27330007254287) q[2];
cx q[2],q[1];
rz(-0.0805141106247902) q[1];
sx q[1];
rz(3.31719948549802) q[1];
sx q[1];
rz(10.3962331771772) q[1];
rz(-0.581412732601166) q[3];
sx q[3];
rz(1.79387799103791) q[3];
sx q[3];
rz(8.04318187235996) q[3];
cx q[3],q[2];
rz(2.4711275100708) q[2];
sx q[2];
rz(5.13107100327546) q[2];
sx q[2];
rz(10.4985816240232) q[2];
rz(0.188875749707222) q[3];
sx q[3];
rz(5.67449942429597) q[3];
sx q[3];
rz(9.04235906004115) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-2.29278659820557) q[0];
sx q[0];
rz(4.50689831574494) q[0];
sx q[0];
rz(10.7474657058637) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(-0.496211916208267) q[1];
sx q[1];
rz(5.19521132309968) q[1];
sx q[1];
rz(11.4876947164456) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(-0.652064800262451) q[2];
sx q[2];
rz(3.87631014187867) q[2];
sx q[2];
rz(10.7911088228147) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(-0.460729569196701) q[3];
sx q[3];
rz(4.98678580124909) q[3];
sx q[3];
rz(14.5151562452237) q[3];
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
