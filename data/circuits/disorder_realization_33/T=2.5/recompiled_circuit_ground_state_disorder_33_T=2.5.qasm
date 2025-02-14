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
rz(1.71409523487091) q[0];
sx q[0];
rz(3.70358982880647) q[0];
sx q[0];
rz(9.19376698731586) q[0];
rz(0.245698735117912) q[1];
sx q[1];
rz(2.68727711041505) q[1];
sx q[1];
rz(8.1375577211301) q[1];
cx q[1],q[0];
rz(0.923169314861298) q[0];
sx q[0];
rz(1.84186628659303) q[0];
sx q[0];
rz(12.0646636247556) q[0];
rz(0.809148609638214) q[2];
sx q[2];
rz(6.61630860169465) q[2];
sx q[2];
rz(11.2570357084195) q[2];
cx q[2],q[1];
rz(-3.41792249679565) q[1];
sx q[1];
rz(4.61385312874848) q[1];
sx q[1];
rz(7.68699750899478) q[1];
rz(-3.54439306259155) q[3];
sx q[3];
rz(9.2346397956186) q[3];
sx q[3];
rz(6.54560277461215) q[3];
cx q[3],q[2];
rz(2.73705577850342) q[2];
sx q[2];
rz(3.84872856934602) q[2];
sx q[2];
rz(9.78527379631206) q[2];
rz(-1.12458419799805) q[3];
sx q[3];
rz(5.23465362389619) q[3];
sx q[3];
rz(11.2974744796674) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.26038134098053) q[0];
sx q[0];
rz(3.31718047161634) q[0];
sx q[0];
rz(8.7678962111394) q[0];
rz(5.95026397705078) q[1];
sx q[1];
rz(4.23407092888887) q[1];
sx q[1];
rz(8.48961088656589) q[1];
cx q[1],q[0];
rz(0.0121977087110281) q[0];
sx q[0];
rz(3.25047881354625) q[0];
sx q[0];
rz(11.2783528327863) q[0];
rz(-1.23360097408295) q[2];
sx q[2];
rz(3.41342750390107) q[2];
sx q[2];
rz(8.01268622874423) q[2];
cx q[2],q[1];
rz(-0.0249273721128702) q[1];
sx q[1];
rz(4.75186601479585) q[1];
sx q[1];
rz(11.1426063537519) q[1];
rz(0.208053395152092) q[3];
sx q[3];
rz(5.06080809433992) q[3];
sx q[3];
rz(7.89535377024814) q[3];
cx q[3],q[2];
rz(1.83445811271667) q[2];
sx q[2];
rz(4.76753273804719) q[2];
sx q[2];
rz(8.20837936400577) q[2];
rz(-0.52792102098465) q[3];
sx q[3];
rz(5.24416771729524) q[3];
sx q[3];
rz(8.16978929042026) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.61277484893799) q[0];
sx q[0];
rz(2.31544938881929) q[0];
sx q[0];
rz(10.0097679257314) q[0];
rz(-0.124755725264549) q[1];
sx q[1];
rz(3.73745331366593) q[1];
sx q[1];
rz(11.1174987316053) q[1];
cx q[1],q[0];
rz(2.12162685394287) q[0];
sx q[0];
rz(3.14559976675595) q[0];
sx q[0];
rz(8.69327870606586) q[0];
rz(0.0135530065745115) q[2];
sx q[2];
rz(4.78456619580323) q[2];
sx q[2];
rz(6.65619847773715) q[2];
cx q[2],q[1];
rz(-2.60220122337341) q[1];
sx q[1];
rz(3.96126470168168) q[1];
sx q[1];
rz(13.1745042562406) q[1];
rz(-0.262890338897705) q[3];
sx q[3];
rz(1.24407270749147) q[3];
sx q[3];
rz(12.3677870988767) q[3];
cx q[3],q[2];
rz(1.88317322731018) q[2];
sx q[2];
rz(4.14341142972047) q[2];
sx q[2];
rz(10.8708350419919) q[2];
rz(2.38407325744629) q[3];
sx q[3];
rz(4.70153001149232) q[3];
sx q[3];
rz(11.7269825696866) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.36767530441284) q[0];
sx q[0];
rz(2.21254191000993) q[0];
sx q[0];
rz(10.5124372005384) q[0];
rz(-3.51678800582886) q[1];
sx q[1];
rz(4.26697841485078) q[1];
sx q[1];
rz(10.3849981784742) q[1];
cx q[1],q[0];
rz(-0.692735731601715) q[0];
sx q[0];
rz(1.61052838166291) q[0];
sx q[0];
rz(9.15949687956973) q[0];
rz(1.60069489479065) q[2];
sx q[2];
rz(4.81179931958253) q[2];
sx q[2];
rz(7.9955750465314) q[2];
cx q[2],q[1];
rz(4.61598062515259) q[1];
sx q[1];
rz(4.28397396405274) q[1];
sx q[1];
rz(15.2147340536039) q[1];
rz(1.23635518550873) q[3];
sx q[3];
rz(5.9986573775583) q[3];
sx q[3];
rz(13.2207507848661) q[3];
cx q[3],q[2];
rz(3.34756922721863) q[2];
sx q[2];
rz(1.87418595154817) q[2];
sx q[2];
rz(7.92721471785709) q[2];
rz(2.28026723861694) q[3];
sx q[3];
rz(3.84437433083589) q[3];
sx q[3];
rz(6.74027345179721) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(3.30615925788879) q[0];
sx q[0];
rz(5.61213866074617) q[0];
sx q[0];
rz(8.82049039601489) q[0];
rz(1.14456486701965) q[1];
sx q[1];
rz(4.62760952313478) q[1];
sx q[1];
rz(12.9888357877652) q[1];
cx q[1],q[0];
rz(-1.23289322853088) q[0];
sx q[0];
rz(2.95412033994729) q[0];
sx q[0];
rz(5.83195707797214) q[0];
rz(0.572117149829865) q[2];
sx q[2];
rz(5.2002102454477) q[2];
sx q[2];
rz(10.7830327510755) q[2];
cx q[2],q[1];
rz(4.62700128555298) q[1];
sx q[1];
rz(2.25861689646775) q[1];
sx q[1];
rz(9.33138510435029) q[1];
rz(-2.29759216308594) q[3];
sx q[3];
rz(5.79114738305146) q[3];
sx q[3];
rz(12.1316601991574) q[3];
cx q[3],q[2];
rz(-1.82332384586334) q[2];
sx q[2];
rz(0.650899799662181) q[2];
sx q[2];
rz(8.76852676867648) q[2];
rz(4.67240619659424) q[3];
sx q[3];
rz(4.41143396695191) q[3];
sx q[3];
rz(8.84404180049106) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.736064732074738) q[0];
sx q[0];
rz(5.27141395409638) q[0];
sx q[0];
rz(11.941390490524) q[0];
rz(2.38963532447815) q[1];
sx q[1];
rz(1.06869211991365) q[1];
sx q[1];
rz(11.0312853813092) q[1];
cx q[1],q[0];
rz(2.29217553138733) q[0];
sx q[0];
rz(4.42348793347413) q[0];
sx q[0];
rz(12.6177158117215) q[0];
rz(0.89959329366684) q[2];
sx q[2];
rz(3.97280648549134) q[2];
sx q[2];
rz(12.3569235563199) q[2];
cx q[2],q[1];
rz(1.24811637401581) q[1];
sx q[1];
rz(-1.79915889899199) q[1];
sx q[1];
rz(8.48345855473682) q[1];
rz(-0.725562632083893) q[3];
sx q[3];
rz(2.40143254597718) q[3];
sx q[3];
rz(12.1021547079007) q[3];
cx q[3],q[2];
rz(2.95251655578613) q[2];
sx q[2];
rz(4.92529729207093) q[2];
sx q[2];
rz(11.6371965169828) q[2];
rz(0.825169861316681) q[3];
sx q[3];
rz(1.6301830132776) q[3];
sx q[3];
rz(11.4639582395475) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-2.58509421348572) q[0];
sx q[0];
rz(0.806990536051341) q[0];
sx q[0];
rz(11.1919164419095) q[0];
rz(-2.63116860389709) q[1];
sx q[1];
rz(8.63428846200044) q[1];
sx q[1];
rz(10.047987496845) q[1];
cx q[1],q[0];
rz(1.57308149337769) q[0];
sx q[0];
rz(3.50911173422868) q[0];
sx q[0];
rz(12.9997985124509) q[0];
rz(6.47896671295166) q[2];
sx q[2];
rz(4.21136120160157) q[2];
sx q[2];
rz(9.93403265475436) q[2];
cx q[2],q[1];
rz(3.67314338684082) q[1];
sx q[1];
rz(1.71677127678926) q[1];
sx q[1];
rz(14.0246047735135) q[1];
rz(-1.78904378414154) q[3];
sx q[3];
rz(8.66520133812959) q[3];
sx q[3];
rz(10.8572345733564) q[3];
cx q[3],q[2];
rz(1.09065425395966) q[2];
sx q[2];
rz(3.9568294604593) q[2];
sx q[2];
rz(10.5920781850736) q[2];
rz(0.0226315706968307) q[3];
sx q[3];
rz(5.76206436951692) q[3];
sx q[3];
rz(11.437443470947) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.24564328789711) q[0];
sx q[0];
rz(4.27679076989228) q[0];
sx q[0];
rz(7.30050680636569) q[0];
rz(-1.74748778343201) q[1];
sx q[1];
rz(3.35269028146798) q[1];
sx q[1];
rz(10.0523273110311) q[1];
cx q[1],q[0];
rz(0.16495569050312) q[0];
sx q[0];
rz(2.40584877331788) q[0];
sx q[0];
rz(11.9541764020841) q[0];
rz(3.05701017379761) q[2];
sx q[2];
rz(5.38990965683992) q[2];
sx q[2];
rz(7.81123051642581) q[2];
cx q[2],q[1];
rz(1.7146475315094) q[1];
sx q[1];
rz(4.94996658165986) q[1];
sx q[1];
rz(12.7680785417478) q[1];
rz(0.340235620737076) q[3];
sx q[3];
rz(1.12005320389802) q[3];
sx q[3];
rz(11.8451981305997) q[3];
cx q[3],q[2];
rz(2.51581859588623) q[2];
sx q[2];
rz(2.50675824482972) q[2];
sx q[2];
rz(9.8355348765771) q[2];
rz(3.09138941764832) q[3];
sx q[3];
rz(7.53610483010346) q[3];
sx q[3];
rz(9.34936617910072) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.55897665023804) q[0];
sx q[0];
rz(5.3863309939676) q[0];
sx q[0];
rz(12.2974531412046) q[0];
rz(-0.310026317834854) q[1];
sx q[1];
rz(4.22863689263398) q[1];
sx q[1];
rz(9.55478776096507) q[1];
cx q[1],q[0];
rz(-0.631424725055695) q[0];
sx q[0];
rz(5.03485170205171) q[0];
sx q[0];
rz(10.7624767780225) q[0];
rz(1.92218458652496) q[2];
sx q[2];
rz(1.23415020306642) q[2];
sx q[2];
rz(17.0075878858487) q[2];
cx q[2],q[1];
rz(4.39062452316284) q[1];
sx q[1];
rz(4.55140498478944) q[1];
sx q[1];
rz(10.4294646739881) q[1];
rz(0.678094506263733) q[3];
sx q[3];
rz(1.16962757905061) q[3];
sx q[3];
rz(13.217792248718) q[3];
cx q[3],q[2];
rz(-2.20760631561279) q[2];
sx q[2];
rz(3.90834847291047) q[2];
sx q[2];
rz(15.8045367956082) q[2];
rz(4.64542579650879) q[3];
sx q[3];
rz(4.47890988190705) q[3];
sx q[3];
rz(8.62507817744418) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.128295168280602) q[0];
sx q[0];
rz(2.95335628290708) q[0];
sx q[0];
rz(13.0994038343351) q[0];
rz(3.17786526679993) q[1];
sx q[1];
rz(0.782500418024608) q[1];
sx q[1];
rz(14.0184082746427) q[1];
cx q[1],q[0];
rz(3.32140946388245) q[0];
sx q[0];
rz(3.84748199780519) q[0];
sx q[0];
rz(11.962580895416) q[0];
rz(3.03259491920471) q[2];
sx q[2];
rz(5.02985611756379) q[2];
sx q[2];
rz(6.26550433634921) q[2];
cx q[2],q[1];
rz(-4.37987184524536) q[1];
sx q[1];
rz(3.5022408088022) q[1];
sx q[1];
rz(10.9106747865598) q[1];
rz(1.39489364624023) q[3];
sx q[3];
rz(5.81131187279756) q[3];
sx q[3];
rz(9.05136713980838) q[3];
cx q[3],q[2];
rz(-2.42862462997437) q[2];
sx q[2];
rz(2.40981009800965) q[2];
sx q[2];
rz(12.599002814285) q[2];
rz(-2.29643583297729) q[3];
sx q[3];
rz(5.05655017693574) q[3];
sx q[3];
rz(10.5049844741742) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.478584825992584) q[0];
sx q[0];
rz(1.75967315037782) q[0];
sx q[0];
rz(11.6238672494809) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(2.32803440093994) q[1];
sx q[1];
rz(4.76537254651124) q[1];
sx q[1];
rz(8.78049484490558) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(1.5004791021347) q[2];
sx q[2];
rz(4.33735290368135) q[2];
sx q[2];
rz(8.86563143729373) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(1.98704385757446) q[3];
sx q[3];
rz(1.83768108685548) q[3];
sx q[3];
rz(11.7202703714292) q[3];
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
