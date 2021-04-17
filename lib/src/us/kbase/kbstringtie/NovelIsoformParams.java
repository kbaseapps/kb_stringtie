
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
 * <p>Original spec-file type: NovelIsoformParams</p>
 * <pre>
 * stringtie_genome_name: name for the new genome including novel transcripts
 * transcript_label: prefix for the name of the output transcripts
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "label",
    "stringtie_genome_name"
})
public class NovelIsoformParams {

    @JsonProperty("label")
    private String label;
    @JsonProperty("stringtie_genome_name")
    private String stringtieGenomeName;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("label")
    public String getLabel() {
        return label;
    }

    @JsonProperty("label")
    public void setLabel(String label) {
        this.label = label;
    }

    public NovelIsoformParams withLabel(String label) {
        this.label = label;
        return this;
    }

    @JsonProperty("stringtie_genome_name")
    public String getStringtieGenomeName() {
        return stringtieGenomeName;
    }

    @JsonProperty("stringtie_genome_name")
    public void setStringtieGenomeName(String stringtieGenomeName) {
        this.stringtieGenomeName = stringtieGenomeName;
    }

    public NovelIsoformParams withStringtieGenomeName(String stringtieGenomeName) {
        this.stringtieGenomeName = stringtieGenomeName;
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
        return ((((((("NovelIsoformParams"+" [label=")+ label)+", stringtieGenomeName=")+ stringtieGenomeName)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
