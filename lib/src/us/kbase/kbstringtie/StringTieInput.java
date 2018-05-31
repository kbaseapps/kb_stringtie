
package us.kbase.kbstringtie;

import java.util.HashMap;
import java.util.Map;
import javax.annotation.Generated;
import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;


/**
 * <p>Original spec-file type: StringTieInput</p>
 * <pre>
 * required params:
 * alignment_object_ref: Alignment or AlignmentSet object reference
 * workspace_name: the name of the workspace it gets saved to
 * expression_set_suffix: suffix append to expression set object name
 * expression_suffix: suffix append to expression object name
 * mode: one of ['normal', 'merge', 'novel_isoform']
 * optional params:
 * num_threads: number of processing threads
 * junction_base: junctions that don't have spliced reads
 * junction_coverage: junction coverage
 * disable_trimming: disables trimming at the ends of the assembled transcripts
 * min_locus_gap_sep_value: minimum locus gap separation value
 * ballgown_mode: enables the output of Ballgown input table files
 * skip_reads_with_no_ref: reads with no reference will be skipped
 * novel_isoforms: output expression matrices with novel isoforms
 * maximum_fraction: maximum fraction of muliple-location-mapped reads
 * min_length: minimum length allowed for the predicted transcripts
 * min_read_coverage: minimum input transcript coverage
 * min_isoform_abundance: minimum isoform abundance
 * ref: http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "alignment_object_ref",
    "workspace_name",
    "expression_set_suffix",
    "expression_suffix",
    "num_threads",
    "junction_base",
    "junction_coverage",
    "disable_trimming",
    "min_locus_gap_sep_value",
    "ballgown_mode",
    "skip_reads_with_no_ref",
    "novel_isoforms",
    "maximum_fraction",
    "min_length",
    "min_read_coverage",
    "min_isoform_abundance"
})
public class StringTieInput {

    @JsonProperty("alignment_object_ref")
    private String alignmentObjectRef;
    @JsonProperty("workspace_name")
    private String workspaceName;
    @JsonProperty("expression_set_suffix")
    private String expressionSetSuffix;
    @JsonProperty("expression_suffix")
    private String expressionSuffix;
    @JsonProperty("num_threads")
    private Long numThreads;
    @JsonProperty("junction_base")
    private Long junctionBase;
    @JsonProperty("junction_coverage")
    private Double junctionCoverage;
    @JsonProperty("disable_trimming")
    private Long disableTrimming;
    @JsonProperty("min_locus_gap_sep_value")
    private Long minLocusGapSepValue;
    @JsonProperty("ballgown_mode")
    private Long ballgownMode;
    @JsonProperty("skip_reads_with_no_ref")
    private Long skipReadsWithNoRef;
    /**
     * <p>Original spec-file type: NovelIsoformParams</p>
     * <pre>
     * stringtie_genome_name: name for the new genome including novel transcripts
     * transcript_label: prefix for the name of the output transcripts
     * </pre>
     * 
     */
    @JsonProperty("novel_isoforms")
    private NovelIsoformParams novelIsoforms;
    @JsonProperty("maximum_fraction")
    private Double maximumFraction;
    @JsonProperty("min_length")
    private Long minLength;
    @JsonProperty("min_read_coverage")
    private Double minReadCoverage;
    @JsonProperty("min_isoform_abundance")
    private Double minIsoformAbundance;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("alignment_object_ref")
    public String getAlignmentObjectRef() {
        return alignmentObjectRef;
    }

    @JsonProperty("alignment_object_ref")
    public void setAlignmentObjectRef(String alignmentObjectRef) {
        this.alignmentObjectRef = alignmentObjectRef;
    }

    public StringTieInput withAlignmentObjectRef(String alignmentObjectRef) {
        this.alignmentObjectRef = alignmentObjectRef;
        return this;
    }

    @JsonProperty("workspace_name")
    public String getWorkspaceName() {
        return workspaceName;
    }

    @JsonProperty("workspace_name")
    public void setWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
    }

    public StringTieInput withWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
        return this;
    }

    @JsonProperty("expression_set_suffix")
    public String getExpressionSetSuffix() {
        return expressionSetSuffix;
    }

    @JsonProperty("expression_set_suffix")
    public void setExpressionSetSuffix(String expressionSetSuffix) {
        this.expressionSetSuffix = expressionSetSuffix;
    }

    public StringTieInput withExpressionSetSuffix(String expressionSetSuffix) {
        this.expressionSetSuffix = expressionSetSuffix;
        return this;
    }

    @JsonProperty("expression_suffix")
    public String getExpressionSuffix() {
        return expressionSuffix;
    }

    @JsonProperty("expression_suffix")
    public void setExpressionSuffix(String expressionSuffix) {
        this.expressionSuffix = expressionSuffix;
    }

    public StringTieInput withExpressionSuffix(String expressionSuffix) {
        this.expressionSuffix = expressionSuffix;
        return this;
    }

    @JsonProperty("num_threads")
    public Long getNumThreads() {
        return numThreads;
    }

    @JsonProperty("num_threads")
    public void setNumThreads(Long numThreads) {
        this.numThreads = numThreads;
    }

    public StringTieInput withNumThreads(Long numThreads) {
        this.numThreads = numThreads;
        return this;
    }

    @JsonProperty("junction_base")
    public Long getJunctionBase() {
        return junctionBase;
    }

    @JsonProperty("junction_base")
    public void setJunctionBase(Long junctionBase) {
        this.junctionBase = junctionBase;
    }

    public StringTieInput withJunctionBase(Long junctionBase) {
        this.junctionBase = junctionBase;
        return this;
    }

    @JsonProperty("junction_coverage")
    public Double getJunctionCoverage() {
        return junctionCoverage;
    }

    @JsonProperty("junction_coverage")
    public void setJunctionCoverage(Double junctionCoverage) {
        this.junctionCoverage = junctionCoverage;
    }

    public StringTieInput withJunctionCoverage(Double junctionCoverage) {
        this.junctionCoverage = junctionCoverage;
        return this;
    }

    @JsonProperty("disable_trimming")
    public Long getDisableTrimming() {
        return disableTrimming;
    }

    @JsonProperty("disable_trimming")
    public void setDisableTrimming(Long disableTrimming) {
        this.disableTrimming = disableTrimming;
    }

    public StringTieInput withDisableTrimming(Long disableTrimming) {
        this.disableTrimming = disableTrimming;
        return this;
    }

    @JsonProperty("min_locus_gap_sep_value")
    public Long getMinLocusGapSepValue() {
        return minLocusGapSepValue;
    }

    @JsonProperty("min_locus_gap_sep_value")
    public void setMinLocusGapSepValue(Long minLocusGapSepValue) {
        this.minLocusGapSepValue = minLocusGapSepValue;
    }

    public StringTieInput withMinLocusGapSepValue(Long minLocusGapSepValue) {
        this.minLocusGapSepValue = minLocusGapSepValue;
        return this;
    }

    @JsonProperty("ballgown_mode")
    public Long getBallgownMode() {
        return ballgownMode;
    }

    @JsonProperty("ballgown_mode")
    public void setBallgownMode(Long ballgownMode) {
        this.ballgownMode = ballgownMode;
    }

    public StringTieInput withBallgownMode(Long ballgownMode) {
        this.ballgownMode = ballgownMode;
        return this;
    }

    @JsonProperty("skip_reads_with_no_ref")
    public Long getSkipReadsWithNoRef() {
        return skipReadsWithNoRef;
    }

    @JsonProperty("skip_reads_with_no_ref")
    public void setSkipReadsWithNoRef(Long skipReadsWithNoRef) {
        this.skipReadsWithNoRef = skipReadsWithNoRef;
    }

    public StringTieInput withSkipReadsWithNoRef(Long skipReadsWithNoRef) {
        this.skipReadsWithNoRef = skipReadsWithNoRef;
        return this;
    }

    /**
     * <p>Original spec-file type: NovelIsoformParams</p>
     * <pre>
     * stringtie_genome_name: name for the new genome including novel transcripts
     * transcript_label: prefix for the name of the output transcripts
     * </pre>
     * 
     */
    @JsonProperty("novel_isoforms")
    public NovelIsoformParams getNovelIsoforms() {
        return novelIsoforms;
    }

    /**
     * <p>Original spec-file type: NovelIsoformParams</p>
     * <pre>
     * stringtie_genome_name: name for the new genome including novel transcripts
     * transcript_label: prefix for the name of the output transcripts
     * </pre>
     * 
     */
    @JsonProperty("novel_isoforms")
    public void setNovelIsoforms(NovelIsoformParams novelIsoforms) {
        this.novelIsoforms = novelIsoforms;
    }

    public StringTieInput withNovelIsoforms(NovelIsoformParams novelIsoforms) {
        this.novelIsoforms = novelIsoforms;
        return this;
    }

    @JsonProperty("maximum_fraction")
    public Double getMaximumFraction() {
        return maximumFraction;
    }

    @JsonProperty("maximum_fraction")
    public void setMaximumFraction(Double maximumFraction) {
        this.maximumFraction = maximumFraction;
    }

    public StringTieInput withMaximumFraction(Double maximumFraction) {
        this.maximumFraction = maximumFraction;
        return this;
    }

    @JsonProperty("min_length")
    public Long getMinLength() {
        return minLength;
    }

    @JsonProperty("min_length")
    public void setMinLength(Long minLength) {
        this.minLength = minLength;
    }

    public StringTieInput withMinLength(Long minLength) {
        this.minLength = minLength;
        return this;
    }

    @JsonProperty("min_read_coverage")
    public Double getMinReadCoverage() {
        return minReadCoverage;
    }

    @JsonProperty("min_read_coverage")
    public void setMinReadCoverage(Double minReadCoverage) {
        this.minReadCoverage = minReadCoverage;
    }

    public StringTieInput withMinReadCoverage(Double minReadCoverage) {
        this.minReadCoverage = minReadCoverage;
        return this;
    }

    @JsonProperty("min_isoform_abundance")
    public Double getMinIsoformAbundance() {
        return minIsoformAbundance;
    }

    @JsonProperty("min_isoform_abundance")
    public void setMinIsoformAbundance(Double minIsoformAbundance) {
        this.minIsoformAbundance = minIsoformAbundance;
    }

    public StringTieInput withMinIsoformAbundance(Double minIsoformAbundance) {
        this.minIsoformAbundance = minIsoformAbundance;
        return this;
    }

    @JsonAnyGetter
    public Map<String, Object> getAdditionalProperties() {
        return this.additionalProperties;
    }

    @JsonAnySetter
    public void setAdditionalProperties(String name, Object value) {
        this.additionalProperties.put(name, value);
    }

    @Override
    public String toString() {
        return ((((((((((((((((((((((((((((((((((("StringTieInput"+" [alignmentObjectRef=")+ alignmentObjectRef)+", workspaceName=")+ workspaceName)+", expressionSetSuffix=")+ expressionSetSuffix)+", expressionSuffix=")+ expressionSuffix)+", numThreads=")+ numThreads)+", junctionBase=")+ junctionBase)+", junctionCoverage=")+ junctionCoverage)+", disableTrimming=")+ disableTrimming)+", minLocusGapSepValue=")+ minLocusGapSepValue)+", ballgownMode=")+ ballgownMode)+", skipReadsWithNoRef=")+ skipReadsWithNoRef)+", novelIsoforms=")+ novelIsoforms)+", maximumFraction=")+ maximumFraction)+", minLength=")+ minLength)+", minReadCoverage=")+ minReadCoverage)+", minIsoformAbundance=")+ minIsoformAbundance)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
